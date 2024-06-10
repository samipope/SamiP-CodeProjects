import java.io.IOException;
import java.net.ServerSocket;
import java.net.Socket;

public class Server {

    private final ServerSocket serverSocket;

    public Server(ServerSocket ss) {
        serverSocket = ss;
    }

    public void runServer() {
        Socket socket;
        while (true) {
            try {
                socket = serverSocket.accept();
                Thread thread = new Thread(new ConnectionHandler(socket));
                thread.start();
            } catch (IOException e) {
                System.out.println(e.getMessage());
                e.printStackTrace();
            }
        }
    }
}