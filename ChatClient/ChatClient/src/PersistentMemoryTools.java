import java.io.*;
import java.util.ArrayList;
import java.util.Map;
import java.util.Scanner;

// read history messages of certain room from local files
// write new messages to local files to retain history messages
public class PersistentMemoryTools {

    private static String getMemoryFileName(String room)
    {
        return "server/memory/" + room + ".txt";
    }

    private static boolean isRoomRecorded(String room)
    {
        File file = new File(getMemoryFileName(room));
        boolean isRecorded = true;
        try {
            new FileInputStream(file);
        } catch (FileNotFoundException e) {
            isRecorded = false;
        }
        return isRecorded;
    }

    private static void createMemoryFileForRoom(String room) throws IOException {
        File directory = new File("server/memory");
        if (!directory.exists()) {
            directory.mkdirs(); // This method creates the directory along with all necessary parent directories.
        }
        File roomFile = new File(directory, room + ".txt");
        try {
            if (!roomFile.exists()) {
                roomFile.createNewFile(); // This creates the file if it does not exist.
            }
            // Now you can safely write to the file
            FileWriter fw = new FileWriter(roomFile, true); // 'true' to append to the file
            // Use FileWriter to write to the file...
        } catch (IOException e) {
            e.printStackTrace();
            throw new RuntimeException(e);
        }
    }

    private static ArrayList<String> readMemoryFileForRoom(String room) throws FileNotFoundException {
        File file = new File(getMemoryFileName(room));
        FileInputStream fileInputStream = new FileInputStream(file);
        Scanner scanner = new Scanner(fileInputStream);
        ArrayList<String> memory = new ArrayList<>();
        while (scanner.hasNextLine())
        {
            String line = scanner.nextLine();
            memory.add(line);
        }
        return memory;
    }

    public static ArrayList<String> getMessageHistoryOfRoom(String room) throws IOException {
        boolean isRoomRecorded = isRoomRecorded(room);
        if(!isRoomRecorded)
        {
            createMemoryFileForRoom(room);
        }
        return readMemoryFileForRoom(room);
    }

    public static void addMessageToMemoryFile(Map<String, String> json) throws IOException {
        String room = json.get("room");
        String message = WebSocketTools.stringifyJSON(json) + "\n";
        boolean isRoomRecorded = isRoomRecorded(room);
        if(!isRoomRecorded)
        {
            createMemoryFileForRoom(room);
        }
        File file = new File(getMemoryFileName(room));
        FileWriter fileWriter = new FileWriter(file, true);
        fileWriter.write(message);
        fileWriter.close();
    }
}