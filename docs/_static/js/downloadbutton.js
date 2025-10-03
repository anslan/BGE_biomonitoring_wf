document.addEventListener("DOMContentLoaded", () => {
  document.querySelectorAll("div.highlight").forEach((block) => {
    let pre = block.querySelector("pre");
    if (!pre) return;

    // Create button
    let button = document.createElement("button");
    button.className = "downloadbtn";
    button.title = "Download code";

    // Unicode  (you can swap for an SVG icon if you want)
    button.innerHTML = "?";

    // Add click handler
    button.addEventListener("click", () => {
      let code = pre.innerText;
      let blob = new Blob([code], { type: "text/plain" });
      let link = document.createElement("a");
      link.href = URL.createObjectURL(blob);
      link.download = "script.txt"; // you could vary this by language, etc.
      link.click();
    });

    // Insert button into the block
    block.style.position = "relative";
    block.appendChild(button);
  });
});
