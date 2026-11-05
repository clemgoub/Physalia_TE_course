# Physalia TE course

Handouts for the practicals of the Physalia course on transposable elements.

| Directory | Practical |
|---|---|
| [`Practical3/`](Practical3/) | Genotyping of mobile element insertions (MEIs) in a TE graph-genome with GraffiTE 1.1 |
| [`Practical4/`](Practical4/) | Differential expression of TEs with TEtranscripts |
| [`render/`](render/) | Script that builds the HTML handouts from the Markdown |

## Handouts to open in a browser

- Practical 3, student handout: https://clemgoub.github.io/Physalia_TE_course/Practical3/Practical3_GraffiTE_v1.1.html
- Practical 3, instructor run sheet: https://clemgoub.github.io/Physalia_TE_course/Practical3/INSTRUCTOR_practical3.html
- Practical 4, student handout: https://clemgoub.github.io/Physalia_TE_course/Practical4/Practical4_TEtranscripts.html
- Practical 4, instructor run sheet: https://clemgoub.github.io/Physalia_TE_course/Practical4/INSTRUCTOR_practical4.html

## The course server

The commands in the handouts run on the course server, where the software and the data are
installed (`/opt/graffite`, `/data/practical3_graffite`, `/data/practical4_tetranscripts`). The server
address and the key to log in are not in the handouts. The instructors give them to the students on the day of the class.

## Building the HTML handouts

The HTML files in this repository were built with [pandoc](https://pandoc.org) 3.7:

```bash
render/render.sh Practical3/Practical3_GraffiTE_v1.1.md "Practical 3: GraffiTE" Practical3/Practical3_GraffiTE_v1.1.html
```

The HTML file embeds the figures, the images and the Mermaid library, so it opens without an
internet connection.
