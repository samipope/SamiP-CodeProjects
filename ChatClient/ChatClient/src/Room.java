import java.io.IOException;
import java.net.Socket;
import java.util.*;

// manage global room states
// track active users & sockets of certain room
public class Room {
    private final String roomName;
    private final Set<String> activeUsers;
    private final Set<Socket> activeSockets;
    private static final Set<Room> roomSet = new HashSet<>();
    private final ArrayList<Map<String, String>> messageQueue;

    Room(String roomName) throws IOException {
        this.roomName = roomName;
        activeUsers = new HashSet<>();
        activeSockets = new HashSet<>();
        messageQueue = new ArrayList<>();
        readMessageQueueFromMemory();
    }

    private void readMessageQueueFromMemory() throws IOException {
        ArrayList<String> messageStrings = PersistentMemoryTools.getMessageHistoryOfRoom(roomName);
        for (String messageString: messageStrings) {
            Map<String, String> json = WebSocketTools.parseJSON(messageString);
            messageQueue.add(json);
        }
    }


    public String getRoomName() {
        return roomName;
    }

    private static synchronized boolean hasRoom(String roomName) {
        for (Room room : roomSet) {
            String currentRoomName = room.getRoomName();
            if (currentRoomName.equals(roomName)) {
                return true;
            }
        }
        return false;
    }

    public static synchronized Room getRoom(String roomName) throws IOException {
        if(hasRoom(roomName))
        {
            for (Room room : roomSet) {
                if (room.getRoomName().equals(roomName)) {
                    return room;
                }
            }
        }
        System.out.println("add room : " + roomName);
        Room newRoom = new Room(roomName);
        roomSet.add(newRoom);
        return newRoom;
    }

    public synchronized void addUser(String user)
    {
        if(activeUsers.contains(user))
        {
            throw new RuntimeException("User already joined the room!");
        }
        System.out.println("add user: " + user);
        activeUsers.add(user);
    }

    public synchronized void removeUser(String user)
    {
        System.out.println("remove user: " + user);
        activeUsers.remove(user);
    }

    public synchronized void addSocket(Socket socket)
    {
        System.out.println("add socket: " + socket);
        activeSockets.add(socket);
    }

    public synchronized void removeSocket(Socket socket)
    {
        System.out.println("remove socket: " + socket);
        activeSockets.remove(socket);
    }

    public synchronized Set<Socket> getActiveSockets()
    {
        return activeSockets;
    }

    private boolean isUserInRoom(String user)
    {
        return activeUsers.contains(user);
    }

    public static synchronized Room getRoomByUser(String user)
    {
        for (Room room: roomSet) {
            if(room.isUserInRoom(user))
            {
                return room;
            }
        }
        throw new RuntimeException("User " + user +" is not in any room!");
    }

    public synchronized void addMessage(Map<String, String> message) throws IOException {
        System.out.println("add message: " + WebSocketTools.stringifyJSON(message));
        messageQueue.add(message);
        PersistentMemoryTools.addMessageToMemoryFile(message);
    }

    public synchronized ArrayList<Map<String, String>> getMessageQueue()
    {
        return messageQueue;
    }

    public boolean canAddUser(String user) {
        return !activeUsers.contains(user);
    }
}