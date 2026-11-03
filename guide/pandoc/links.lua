-- Pandoc Lua filter for the dgretl guide:
--  * links to other markdown files (foo.md, docs/x/y.md#frag) point to the html files instead;
--  * the page title (pagetitle) is the text of the first level-1 heading.

function Link(el)
  local t = el.target
  if t:match("^%a[%w+.-]*:") or t:match("^#") then
    return nil            -- an absolute url or a fragment of this page
  end
  local path, frag = t:match("^([^#]*)(#?.*)$")
  if path:match("%.md$") then
    el.target = path:gsub("%.md$", ".html") .. frag
    return el
  end
  return nil
end

function Pandoc(doc)
  if doc.meta.pagetitle == nil then
    for _, block in ipairs(doc.blocks) do
      if block.t == "Header" and block.level == 1 then
        doc.meta.pagetitle = pandoc.MetaString(pandoc.utils.stringify(block.content))
        break
      end
    end
  end
  return doc
end
