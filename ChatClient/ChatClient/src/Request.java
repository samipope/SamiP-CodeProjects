import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.net.Socket;
import java.nio.charset.StandardCharsets;
import java.util.HashMap;
import java.util.Map;
import java.util.Scanner;

// read & parse incoming requests
public class Request {
    private final Socket socket;
    private final Map<String, String> requestHeader = new HashMap<>();
    private String fileName = "index.html";
    private Boolean isFileValid;
    private String httpVersion = "";

    private boolean isWebSocket = false;

    private String websocketRequestKey = "";

    public Request(Socket s) {
        socket = s;
    }

    public String getFilePath(String fn) {
      //  String basePath = "../../resources";
       String basePath = "/Users/samanthapope/PortfolioGithub/ChatClient/ChatClient/resources/";
       // String basePath = "/Users/samanthapope/MSD/Github/CS6011/Day20/ChatClient/resources";
        return basePath + fn; // Ensure 'fn' begins with '/' to construct the path correctly.

    }

    private void checkIfFileExists() {
        File file = new File(getFilePath(fileName));
        isFileValid = file.exists() && !file.isDirectory();
    }

    private void printKeyValueMap() {
        requestHeader.forEach((key, value) -> System.out.println(key + " : " + value));
    }

    private void readRequestHeader(Scanner requestHeaderScanner) {
        // read the first line for the fileName
        // e.g. GET /tutorials/other/top-20-mysql-best-practices/ HTTP/1.1
        String GETInfo = requestHeaderScanner.nextLine();
        System.out.println(GETInfo + " thread: " + Thread.currentThread().threadId());
        String[] GETInfoArray = GETInfo.split("\s");
        fileName = GETInfoArray[1];
        checkIfFileExists();
        httpVersion = GETInfoArray[2];
        System.out.println(GETInfo);
        // keep reading the following key-value pairs and store them in requestHeader
        while (requestHeaderScanner.hasNextLine()) {
            // e.g. Host: localhost
            String keyValueString = requestHeaderScanner.nextLine();
            // end at the blank line
            if (keyValueString.isEmpty()) {
                break;
            }
            String[] keyValuePair = keyValueString.split(": ");
            String key = keyValuePair[0];
            String value = keyValuePair[1];
            requestHeader.put(key, value);
        }
        isWebSocket = WebSocketTools.isWebSocket(requestHeader);
        if (isWebSocket) {
            websocketRequestKey = requestHeader.get("Sec-WebSocket-Key");
        }
        printKeyValueMap();
    }

    public Boolean getIsFileValid() {
        return isFileValid;
    }

    public String getFileName() {
        return fileName;
    }

    public String getHttpVersion() {
        return httpVersion;
    }

    public void getRequest() {
        InputStream requestStream = null;
        try {
            requestStream = socket.getInputStream();
        } catch (IOException e) {
            System.out.println(e.getMessage());
            e.printStackTrace();
        }
        assert requestStream != null;
        Scanner requestHeaderScanner = new Scanner(requestStream, StandardCharsets.UTF_8);
        if (requestHeaderScanner.hasNextLine()) {
            readRequestHeader(requestHeaderScanner);
        }
    }


    public boolean getIsWebSocket() {
        return isWebSocket;
    }

    public String getWebsocketRequestKey() {
        if (!isWebSocket) {
            throw new RuntimeException("This is not a web socket!");
        }
        return websocketRequestKey;
    }
}