import java.io.*;
import java.net.Socket;
import java.nio.charset.StandardCharsets;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.*;

// read websocket requests & send websocket responses
public class WebSocketTools {
    private static final String magicString = "258EAFA5-E914-47DA-95CA-C5AB0DC85B11";

    public static String getWebSocketResponseKey(String requestKey) throws NoSuchAlgorithmException {
        String concatenatedString = requestKey + magicString;
        MessageDigest sha1 = MessageDigest.getInstance("SHA-1");
        byte[] sha1Results = sha1.digest(concatenatedString.getBytes(StandardCharsets.UTF_8));
        return Base64.getEncoder().encodeToString(sha1Results);
    }


    public static boolean isWebSocket(Map<String, String> requestHeader) {
        return requestHeader.containsKey("Upgrade") && requestHeader.get("Upgrade").equals("websocket") && requestHeader.containsKey("Sec-WebSocket-Key");
    }

    private static boolean isMasked(ArrayList<Byte> message) {
        return (message.get(1) & 0x80) != 0;
    }


    public static String getWebSocketResponsePayload(ArrayList<Byte> message, long payloadStart, long payloadLength) {
        boolean isMasked = isMasked(message);
        byte[] payload = new byte[(int) payloadLength];
        for (int i = 0; i < payloadLength; i++) {
            payload[i] = message.get((int) (i + payloadStart));
        }
        if (isMasked) {
            // decode the masked message
            byte[] masks = new byte[4];
            for (int i = 0; i < 4; i++) {
                masks[i] = message.get((int) (i + payloadStart - 4));
            }
            for (int i = 0; i < payloadLength; i++) {
                payload[i] = (byte) (payload[i] ^ masks[i % 4]);
            }
        }
        return new String(payload, StandardCharsets.UTF_8);
    }

    private static ArrayList<Byte> getBytesFromNumber(int number, int length) {
        ArrayList<Byte> results = new ArrayList<>();
        for(int i=0; i<length; i++)
        {
            byte current = (byte) (number & 0xff);
            number = number >> 8;
            results.add(current);
        }
        Collections.reverse(results);
        return results;
    }

    public static byte[] getResponseFrame(String response) {
        byte[] payload = response.getBytes(StandardCharsets.UTF_8);
        int payloadLength = payload.length;
        ArrayList<Byte> responseFrame = new ArrayList<>();
        // FIN / RSV*3 / OPCODE
        byte firstByte = (byte) 0x81;
        responseFrame.add(firstByte);
        if (payloadLength <= 125) {
            byte secondByte = (byte) payloadLength;
            responseFrame.add(secondByte);
        } else if (payloadLength <= 0xffff) {
            byte secondByte = 126;
            responseFrame.add(secondByte);
            ArrayList<Byte> lengthBytes = getBytesFromNumber(payloadLength, 2);
            responseFrame.addAll(lengthBytes);
        } else {
            byte secondByte = 127;
            responseFrame.add(secondByte);
            ArrayList<Byte> lengthBytes = getBytesFromNumber(payloadLength, 8);
            responseFrame.addAll(lengthBytes);
        }
        for (byte p :
                payload) {
            responseFrame.add(p);
        }
        byte[] results = new byte[responseFrame.size()];
        for (int i = 0; i < responseFrame.size(); i++) {
            results[i] = responseFrame.get(i);
        }
        return results;
    }

    public static void makeHandshake(Socket socket, Request request) {
        OutputStream out;
        try {
            out = socket.getOutputStream();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        byte[] response;
        try {
            response = ("HTTP/1.1 101 Switching Protocols\r\n"
                    + "Connection: Upgrade\r\n"
                    + "Upgrade: websocket\r\n"
                    + "Sec-WebSocket-Accept: "
                    + WebSocketTools.getWebSocketResponseKey(request.getWebsocketRequestKey())
                    + "\r\n\r\n").getBytes(StandardCharsets.UTF_8);
        } catch (NoSuchAlgorithmException e) {
            throw new RuntimeException(e);
        }
        try {
            out.write(response);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        System.out.println("handshake with the client");
    }

    public static String getRequest(DataInputStream dataInputStream) throws IOException {
        ArrayList<Byte> incomingBytes = new ArrayList<>();
        // read first 2 bytes
        for (int i = 0; i < 2; i++) {
            incomingBytes.add(dataInputStream.readByte());
        }
        boolean isMasked = isMasked(incomingBytes);
        boolean isClosingFrame = isClosingFrame(incomingBytes);
        int byte1PayloadLength = incomingBytes.get(1) & 0x7f;
        long payloadLength;
        int payloadStart;
        if (byte1PayloadLength < 126) {
            // payload length in byte 1
            payloadLength = byte1PayloadLength;
            payloadStart = 2;
        } else if (byte1PayloadLength == 126) {
            // extended payload length in byte 2 ~ 3
            // read next 2 bytes / 1 short
            payloadLength = dataInputStream.readShort();
            short currentPayloadLength = (short) payloadLength;
            for (int i = 0; i < 2; i++) {
                incomingBytes.add((byte) (currentPayloadLength & 0xff));
                currentPayloadLength = (short) (currentPayloadLength >> 8);
            }
            payloadStart = 4;
        } else {
            // extended payload length in byte 2 ~ 9
            // read next 8 bytes / 1 long
            payloadLength = dataInputStream.readLong();
            long currentPayloadLength = payloadLength;
            for (int i = 0; i < 8; i++) {
                incomingBytes.add((byte) (currentPayloadLength & 0xff));
                currentPayloadLength = (short) (currentPayloadLength >> 8);
            }
            payloadStart = 10;
        }
        if (isMasked) {
            // read next 4 bytes as masks
            for (int i = 0; i < 4; i++) {
                incomingBytes.add(dataInputStream.readByte());
            }
            payloadStart += 4;
        }
        // read next payloadLength bytes for the payload
        for (int i = 0; i < payloadLength; i++) {
            incomingBytes.add(dataInputStream.readByte());
        }
        if (isClosingFrame) {
            return "close";
        }
        return getWebSocketResponsePayload(incomingBytes, payloadStart, payloadLength);
    }

    private static boolean isClosingFrame(ArrayList<Byte> incomingBytes) {
        // OPCODE 8 : CLOSE
        return (incomingBytes.get(0) & 0x0f) == 8;
    }

    public static void sendResponse(String response, Socket socket) throws IOException {
        System.out.println("send response: " + response);
        OutputStream out = null;
        try {
            out = socket.getOutputStream();
        } catch (IOException e) {
            System.out.println(e.getMessage());
            e.printStackTrace();
        }
        assert out != null;
        out.write(getResponseFrame(response));
    }


    private static Map<String, String> getResponseJSON(String request) {
        // Assume request format is "action user room message" for messages
        // and "action user room" for join/leave
        String[] incomingRequests = request.split("\\s+", 4); // Split by space, limit to 4 parts
        Map<String, String> responseJSON = new HashMap<>();
        if (incomingRequests.length < 3) {
            responseJSON.put("error", "Invalid request format.");
            return responseJSON;
        }
        String action = incomingRequests[0];
        switch (action) {
            case "join", "leave" -> {
                if (incomingRequests.length == 3) { // Correct length for join/leave
                    responseJSON.put("type", action);
                    responseJSON.put("user", incomingRequests[1]);
                    responseJSON.put("room", incomingRequests[2]);
                } else {
                    responseJSON.put("error", "Invalid join/leave request format.");
                }
            }
            case "message" -> {
                if (incomingRequests.length == 4) { // Expecting 4 parts for message
                    responseJSON.put("type", action);
                    responseJSON.put("user", incomingRequests[1]);
                    responseJSON.put("room", incomingRequests[2]);
                    responseJSON.put("message", incomingRequests[3]); // The actual message text
                } else {
                    responseJSON.put("error", "Invalid message request format.");
                }
            }
            default -> responseJSON.put("error", "Unknown action.");
        }
        return responseJSON;
    }

    public static String stringifyJSON(Map<String, String> json) {
        StringBuilder result = new StringBuilder();
        int i = 0;
        for (String key : json.keySet()) {
            result.append("\"").append(key).append("\":").append("\"").append(json.get(key)).append("\"");
            i++;
            if (i < json.size()) {
                result.append(",");
            }
        }
        return "{" + result + "}";
    }

    private static void handleJoin(Map<String, String> responseJSON, Socket socket) throws IOException {
        Room room = Room.getRoom(responseJSON.get("room"));
        room.addUser(responseJSON.get("user"));
        room.addSocket(socket);
    }

    private static void handleLeave(Map<String, String> responseJSON, Socket socket) throws IOException {
        Room room = Room.getRoom(responseJSON.get("room"));
        room.removeUser(responseJSON.get("user"));
        room.removeSocket(socket);
    }

    public static void broadcastResponse(Map<String, String> responseJSON, Socket socket) throws IOException {
        String JSONString = stringifyJSON(responseJSON);
        Room room = Room.getRoom(responseJSON.get("room"));
        Set<Socket> clients = room.getActiveSockets();
        System.out.println("active clients length: " + clients.size());
        if (Objects.equals(responseJSON.get("type"), "join")) {
            // send all past messages to the new joined client
            ArrayList<Map<String, String>> messageQueue = room.getMessageQueue();
            for (Map<String, String> pastMessage : messageQueue) {
                sendResponse(WebSocketTools.stringifyJSON(pastMessage), socket);
            }
        }
        for (Socket client : clients) {
            if (!Objects.equals(responseJSON.get("type"), "join") || !client.equals(socket)) {
                sendResponse(JSONString, client);
            }
        }
    }

    public static void handleResponse(Socket socket, String request) throws IOException {
        Map<String, String> responseJSON = getResponseJSON(request);
        String type = responseJSON.get("type");

        if (type == null || responseJSON.containsKey("error")) {
            // Log error, send error response to client, or handle as appropriate
            System.out.println("Error processing request: " + responseJSON.getOrDefault("error", "Unknown error"));
            return; // Stop further processing to avoid NullPointerException
        }
        Room room = Room.getRoom(responseJSON.get("room"));
        if (responseJSON.containsKey("room")) {
            String roomName = responseJSON.get("room");
            room = Room.getRoom(roomName);
            if (room == null) {
                System.out.println("Error: Room '" + roomName + "' does not exist.");
                return;
            }
        }

        switch (type) {
            case "join" -> {
                if (room != null && room.canAddUser(responseJSON.get("user"))) {
                    handleJoin(responseJSON, socket);
                } else {
                    // Optionally, send a response back to the client indicating the user cannot join.
                    System.out.println("Error: User '" + responseJSON.get("user") + "' cannot join the room '" + responseJSON.get("room") + "'.");
                }
            }
            case "leave" -> {
                if (room != null) {
                    handleLeave(responseJSON, socket);
                }
            }
            case "message" -> {
                if (room != null) {
                    room.addMessage(responseJSON);
                }
            }
            default -> System.out.println("Error: Unknown action type '" + type + "'.");
        }
        // after handling join/leave/message, broadcast the response if necessary.
        if (room != null && !"join".equals(type) && !"leave".equals(type)) {
            broadcastResponse(responseJSON, socket);
        }

    }

    private static void handleDuplicateUserError(Map<String, String> responseJSON, Socket socket) throws IOException {
        Map<String, String> json = new HashMap<>();
        json.put("type", "error");
        json.put("error", "There is someone called " + responseJSON.get("user") + " in " + responseJSON.get("room") + "! Please choose a different user name.");
        json.put("timestamp", new String(String.valueOf(System.currentTimeMillis())));
        String jsonString = stringifyJSON(json);
        sendResponse(jsonString, socket);
    }

    public static Map<String, String> parseJSON(String str) {
        Map<String, String> json = new HashMap<>();
        String trimStr = str.substring(1, str.length() - 1);
        String[] keyValuePairs = trimStr.split(",");
        for (String keyValuePair : keyValuePairs) {
            String[] keyAndValue = keyValuePair.split(":");
            String key = keyAndValue[0].substring(1, keyAndValue[0].length() - 1);
            String value = keyAndValue[1].substring(1, keyAndValue[1].length() - 1);
            json.put(key, value);
        }
        return json;
    }
}