# Physalia TE course

Handouts for the practicals of the Physalia course on transposable elements.

| Directory | Practical |
|---|---|
| [`Practical3/`](Practical3/) | Genotyping of mobile element insertions (MEIs) in a TE graph-genome with GraffiTE 1.1 |
| [`render/`](render/) | Script that builds the HTML handouts from the Markdown |

Practical 4 (TEtranscripts) will go in `Practical4/`.

## The course server

The commands in the handouts run on the course server, where the software and the data are
installed (`/opt/graffite`, `/data/practical3_graffite`). The server address and the key to log in
are not in the handouts. The instructors give them to the students on the day of the class.

## Building the HTML handouts

The HTML files in this repository were built with [pandoc](https://pandoc.org) 3.7:

```bash
render/render.sh Practical3/Practical3_GraffiTE_v1.1.md "Practical 3: GraffiTE" Practical3/Practical3_GraffiTE_v1.1.html
```

The HTML file embeds the figures, the images and the Mermaid library, so it opens without an
internet connection.
