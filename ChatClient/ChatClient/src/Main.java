import java.net.ServerSocket;
import java.io.IOException;

public class Main {
    public static void main(String[] args) {
        // listen at port 8080
        try (ServerSocket serverSocket = new ServerSocket(8080)) {
            Server server = new Server(serverSocket);
            server.runServer();
        } catch (IOException e) {
            System.out.println(e.getMessage());
            e.printStackTrace();
        }
    }
}

//TODO add in comments at the top of each of the classes about purpose
//TODO fix the