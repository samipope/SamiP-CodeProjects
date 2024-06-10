# Chat Client

## Languages Used
- Java
- HTML
- CSS
- JavaScript


## Features
- **Real-Time Communication**: Enables instant messaging across different browsers, utilizing WebSocket technology for efficient, live data transfer.
- **Concurrency Management**: Manages multiple user connections simultaneously, showcasing skills in handling concurrency and multithreading in Java.
- **Full-Stack Development**: Incorporates both front-end and back-end technologies including HTML, CSS, JavaScript, and Java, demonstrating a comprehensive approach to web application development.
- **Scalable Architecture**: Designed to support multiple users and rooms without degradation in performance, illustrating principles of scalable software architecture.
- **Network Programming**: Extensive use of network protocols (HTTP, WebSocket) to establish and maintain communications, highlighting skills in network programming.
- **User Experience Design**: Ensures a seamless and responsive user interface on various devices and browsers, emphasizing attention to cross-platform compatibility and user experience.
- **JSON Data Handling**: Utilizes JSON for messaging and data interchange between clients and server, showing proficiency in working with data formats and APIs.


<img width="1399" alt="Screenshot 2024-06-10 at 3 10 14 PM" src="https://github.com/samipope/SamiP-CodeProjects/assets/142822253/cd6146d2-053b-46c5-84c9-f3f376a50bd3">


## Project Summary:
The ChatClient project is a sophisticated real-time communication system that allows users to connect and chat across different browsers within designated "rooms". Utilizing WebSockets, the system facilitates bi-directional communication between the server and clients, beginning with an HTTP request that upgrades to a WebSocket connection. This ensures all messages within a room are broadcasted instantly to all participants, enhancing interactive user experience.

In this application, users can join rooms and send messages, which are managed through a combination of JSON parsing and WebSocket communication. The backend, built with Java, handles concurrent connections efficiently using multithreading, illustrating robust handling of network protocols and concurrency. The frontend leverages HTML, CSS, and JavaScript, demonstrating a full-stack development approach.

While developing this project, I was able to deepen my understanding of network programming, concurrency management, and real-time data handling. Working on this project taught me valuable lessons in managing complex communication protocols and enhancing user interactions through modern web technologies.

## How to Run
1. **Start the Server**:
   - Execute `Main.java` to start the server.

2. **Access the Chat Interface**:
   - Open a web browser and navigate to `http://localhost:8080/chat.html`.

3. **Connecting Multiple Clients**:
   - To add another client, open a new browser window or a different browser (e.g.,   Chrome, Firefox) and access the same URL.
   - You can also connect from different devices on the same network by using the host machine's IP address instead of `localhost`, e.g., `http://192.168.x.x:8080/chat.html`.
