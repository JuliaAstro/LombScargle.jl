using DocumenterCitations

pages = ["Home" => "index.md"]


# put bib here for juliaastro.github.io so that the builder for that site can
# pick it up
bib = CitationBibliography(
    joinpath(@__DIR__, "src", "references.bib");
    style=:authoryear,
)

