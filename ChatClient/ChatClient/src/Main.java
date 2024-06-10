import java.net.ServerSocket;
import java.io.IOException;

public class Main {
    public static void main(String[] args) throws IOException {
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
