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
      let code = "";
      
      // Check if line numbers are in a separate table column (Sphinx table format)
      let table = block.querySelector("table.highlighttable");
      if (table) {
        // Get code from the code cell (not the line numbers column)
        let codeCell = table.querySelector("td.code pre");
        if (codeCell) {
          code = codeCell.innerText;
        } else {
          code = pre.innerText;
        }
      } else {
        // Line numbers are inline - clone the pre and remove line number elements
        let preClone = pre.cloneNode(true);
        
        // Remove all line number spans/elements (Sphinx uses span.linenos for line numbers)
        // This includes: .lineno, .linenos, .linenodiv, and any span with these classes
        preClone.querySelectorAll("span.linenos, span.lineno, .lineno, .linenos, .linenodiv").forEach(el => {
          el.remove();
        });
        
        // Extract text content - this preserves the original code structure and indentation
        // innerText handles whitespace correctly and excludes removed elements
        code = preClone.innerText || preClone.textContent;
      }
      
      // Determine file extension
      let filename = "script.txt";
      let trimmedCode = code.trimStart();
      if (trimmedCode.startsWith("#!/bin/bash") || trimmedCode.startsWith("#!/usr/bin/bash") || trimmedCode.startsWith("#!/bin/sh")) {
        filename = "script.sh";
      } else if (trimmedCode.startsWith("#!/usr/bin/Rscript")) {
        filename = "script.R";
      }
      
      let blob = new Blob([code], { type: "text/plain" });
      let link = document.createElement("a");
      link.href = URL.createObjectURL(blob);
      link.download = filename;
      link.click();
    });

    // Insert button
    block.style.position = "relative";
    block.appendChild(button);
  });
});
