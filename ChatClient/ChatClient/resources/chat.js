const username = document.getElementById("Username");
const chatroom = document.getElementById("chatroom");
const message1 = document.getElementById("message");
const sendButton = document.getElementById("sendButton");
const leaveRoomButton = document.getElementById("leaveRoomButton");
const joinRoomButton = document.getElementById("joinRoomButton");
const ws = new WebSocket("ws://localhost:8080");

let wsOpen = false;
let userJoined = false;

ws.onopen = function () {
    wsOpen = true;
};

function sendMessage() {
    if (wsOpen && userJoined) {
        const messageValue = message1.value.trim();
        const chatroomValue = chatroom.value.trim();
        const usernameValue = username.value.trim();
        // Correctly format the message
        ws.send(`message ${usernameValue} ${chatroomValue} ${messageValue}`);
        console.log('message sent: ' + messageValue);
       // displayMessage(usernameValue, messageValue);
    }
}

function joinRoom() {
    if (wsOpen && !userJoined) {
        const chatroomValue = chatroom.value.trim();
        const usernameValue = username.value.trim();
        // Ensure chatroom name is lowercase
        if (!isChatroomValid(chatroomValue)) {
            alert("Chatroom name should be all lowercase!");
            chatroom.focus();
            return;
        }
        ws.send(`join ${usernameValue} ${chatroomValue}`);
        userJoined = true;
        console.log(usernameValue + ' joined ' + chatroomValue);
        displayJoin(usernameValue);
    }
}

function leaveRoom() {
    if (wsOpen && userJoined) {
        const usernameValue = username.value.trim();
        ws.send(`leave ${usernameValue}`);
        userJoined = false;
        console.log(usernameValue + ' left the chat');
        displayLeave(usernameValue);
    }
}

ws.onmessage = function (event) {
    let msgObj = JSON.parse(event.data);
    if (msgObj.type === 'message') {
        displayMessage(msgObj.user, msgObj.message);
    } else if (msgObj.type === 'join') {
        displayJoin(msgObj.user);
    } else if (msgObj.type === 'leave') {
        displayLeave(msgObj.user);
    }
};

function displayMessage(user, message) {
    let messagePar = document.getElementById("Messages");
    let msgElement = document.createElement("p");
    msgElement.textContent = user + ": " + message;
    messagePar.appendChild(msgElement);
}

function displayJoin(user) {
    let peoplePar = document.getElementById("People");
    let pplElement = document.createElement("div");
    pplElement.textContent = user + " has joined the room";
    peoplePar.appendChild(pplElement);
}

function displayLeave(user) {
    let peoplePar = document.getElementById("People");
    let pplElement = document.createElement("div");
    pplElement.textContent = user + " has left the chat";
    peoplePar.appendChild(pplElement);
}

function isChatroomValid(name) {
    // Check if chatroom name is lowercase
    return name === name.toLowerCase();
}

// Add event listeners for sending messages on "Enter"
message1.addEventListener("keypress", (event) => {
    if (event.key === "Enter") {
        sendMessage();
    }
});

// Attach event listeners to buttons
sendButton.addEventListener("click", sendMessage);
joinRoomButton.addEventListener("click", joinRoom);
leaveRoomButton.addEventListener("click", leaveRoom);
