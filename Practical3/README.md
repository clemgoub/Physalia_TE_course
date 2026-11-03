# Practical 3: GraffiTE

| File | Contents |
|---|---|
| [`Practical3_GraffiTE_v1.1.md`](Practical3_GraffiTE_v1.1.md) | Student handout. GitHub shows it with its figures and diagrams. |
| `Practical3_GraffiTE_v1.1.html` | The same handout as a single HTML file. Download it and open it in a browser. |
| [`INSTRUCTOR_practical3.md`](INSTRUCTOR_practical3.md), `.html` | Instructor run sheet: the commands to type live, and how long each took on the course server on 3 November 2026. |
| [`figures/`](figures/) | The cartoons in the handout (SVG), and `make_figures.py`, the script that draws them. |
| [`diagrams/`](diagrams/) | The Mermaid code of the two workflow diagrams. The handout contains the same code: if you change a diagram, change it in both. |

To redraw the figures and rebuild the HTML, from this directory:

```bash
python3 figures/make_figures.py
../render/render.sh Practical3_GraffiTE_v1.1.md "Practical 3: GraffiTE" Practical3_GraffiTE_v1.1.html
```
