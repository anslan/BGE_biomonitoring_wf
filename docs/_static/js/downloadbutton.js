document.addEventListener("DOMContentLoaded", () => {
  document.querySelectorAll("div.highlight pre").forEach((block) => {
    // Create button
    let button = document.createElement("button");
    button.innerText = "⭳ Download";
    button.classList.add("downloadbtn");

    // Insert button
    block.parentNode.insertBefore(button, block);

    // Add click handler
    button.addEventListener("click", () => {
      let code = block.innerText;
      let blob = new Blob([code], {type: "text/plain"});
      let link = document.createElement("a");
      link.href = URL.createObjectURL(blob);
      link.download = "script.txt";
      link.click();
    });
  });
});
