document.addEventListener("DOMContentLoaded", () => {
  document.querySelectorAll("div.highlight").forEach((block) => {
    let pre = block.querySelector("pre");
    if (!pre) return;

    // Create button
    let button = document.createElement("button");
    button.className = "downloadbtn";
    button.title = "Download code";

    // Use a down-arrow SVG (cleaner than Unicode)
    button.innerHTML = `
      <svg xmlns="http://www.w3.org/2000/svg" width="14" height="14" 
           fill="currentColor" viewBox="0 0 16 16">
        <path d="M.5 9.9a.5.5 0 0 1 .5-.5h4V1.5a.5.5 0 0 1 1 0v7.9h4a.5.5 0 0 1 .4.8l-4.5 5a.5.5 0 0 1-.8 0l-4.5-5a.5.5 0 0 1-.1-.3z"/>
      </svg>
    `;

    // Click handler → download content
    button.addEventListener("click", () => {
      let code = pre.innerText;
      let blob = new Blob([code], { type: "text/plain" });
      let link = document.createElement("a");
      link.href = URL.createObjectURL(blob);
      link.download = "script.py"; // default
      link.click();
    });

    // Insert button
    block.style.position = "relative";
    block.appendChild(button);
  });
});
