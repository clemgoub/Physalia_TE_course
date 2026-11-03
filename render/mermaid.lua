-- Turn ```mermaid code blocks into <pre class="mermaid">, which mermaid.js draws in the browser.
local function escape(s)
  return (s:gsub("&", "&amp;"):gsub("<", "&lt;"):gsub(">", "&gt;"))
end
function CodeBlock(el)
  if el.classes[1] == "mermaid" then
    return pandoc.RawBlock("html", '<pre class="mermaid">' .. escape(el.text) .. "</pre>")
  end
end
