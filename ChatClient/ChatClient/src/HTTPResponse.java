import java.io.*;
import java.net.Socket;
import java.nio.file.Files;

// send a HTTP response
public class HTTPResponse {
    private final Socket socket;
    private final Boolean isFileValid;
    private final String fileName;
    private final String httpVersion;
    private final Request request;
    private final String fallback404PageFileName = "/index.html";

    public HTTPResponse(Socket s, Request request) {
        socket = s;
        this.request = request;
        isFileValid = request.getIsFileValid();
        fileName = request.getFileName();
        httpVersion = request.getHttpVersion();
    }

    private String getFinalFilePath() {
        // if the file does not exist, return 404 fall back page
        return request.getFilePath(isFileValid ? fileName : fallback404PageFileName);
    }

    private String getFallBackPageHTML(String errorMessage) {
        return """
                <!doctype html>
                 <html lang="en">
                 <head>
                     <meta charset="utf-8"/>
                     <link rel="icon" href="favicon.webp"/>
                     <meta name="viewport" content="width=device-width,initial-scale=1"/>
                     <meta name="theme-color" content="#000000"/>
                   
                     <title>404 NOT FOUND</title>
                     <link href="404PageStyle.css" rel="stylesheet">
                 </head>
                 <body>
                 <div id="root">
                  
                         <h2> Hmm... What you want is not here</h2>
                         <h3> Error Message From the Server: </h3>
                         <p>
                          """ +
                errorMessage
                +
                """
                        </p>
                            </div>
                        </div>
                        </body>
                        </html>
                        """;
    }

    private void write404fallbackPage() {
        try {
            new FileInputStream(request.getFilePath(fileName));
        } catch (FileNotFoundException e) {
            try {
                FileWriter fileWriter = new FileWriter(request.getFilePath(fallback404PageFileName));
                fileWriter.write(getFallBackPageHTML(e.getMessage()));
                fileWriter.flush();
                fileWriter.close();
            } catch (IOException ex) {
                System.out.println(ex.getMessage());
                ex.printStackTrace();
            }
        }
    }

    private String getStatusCodeInfo() {
        String statusCodeInfo = isFileValid ? "200 OK" : "404 NOT FOUND";
        return httpVersion + " " + statusCodeInfo;
    }

    private String getContentType() {
        File file = new File(getFinalFilePath());
        String fileContentType = "";
        try {
            fileContentType = Files.probeContentType(file.toPath());
        } catch (IOException e) {
            System.out.println("Failed: unknown content of " + fileName);
            e.printStackTrace();
        }
        return "content-type: " + fileContentType;
    }

    private void sendResponseHeader(PrintWriter printWriter) {
        if (!isFileValid) write404fallbackPage();
        printWriter.println(getStatusCodeInfo());
        printWriter.println(getContentType());
        // add blank line to indicate the start of the requested file content
        printWriter.println("");
    }

    private void sendResponseBody(OutputStream socketOutputStream) {
        try {
            InputStream finalFileStream = new FileInputStream(getFinalFilePath());
            finalFileStream.transferTo(socketOutputStream);
            finalFileStream.close();
        } catch (IOException e) {
            // if the file does not exist then return a fallback page to display server error message
            System.out.println(e.getMessage());
            e.printStackTrace();
        }
    }

    private void sendResponse(Socket socket) {
        try {
            OutputStream socketOutputStream = socket.getOutputStream();
            PrintWriter printWriter = new PrintWriter(socketOutputStream, true);
            sendResponseHeader(printWriter);
            //Pause for 1 second
            Thread.sleep(1000);
            sendResponseBody(socketOutputStream);
            printWriter.close();
            socketOutputStream.close();
        } catch (IOException | InterruptedException e) {
            System.out.println(e.getMessage());
            e.printStackTrace();
        }
    }

    public void handleResponse() {
        sendResponse(socket);
    }

}