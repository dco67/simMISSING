-- encadres.lua — rend les encadrés .app (« À l'application ») et .voir
-- (« Ce qu'on doit voir ») dans les trois formats.
--   HTML / EPUB : classes conservées, mises en forme par styles.css
--   PDF        : environnements tcolorbox définis dans preamble.tex

function Div(el)
  if not FORMAT:match("latex") then return nil end
  local env = nil
  if el.classes:includes("app") then env = "encadreapp"
  elseif el.classes:includes("voir") then env = "encadrevoir" end
  if env == nil then return nil end
  local out = { pandoc.RawBlock("latex", "\\begin{" .. env .. "}") }
  for _, b in ipairs(el.content) do table.insert(out, b) end
  table.insert(out, pandoc.RawBlock("latex", "\\end{" .. env .. "}"))
  return out
end
