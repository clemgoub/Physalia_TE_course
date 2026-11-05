# Practical 4: TEtranscripts

| File | Contents |
|---|---|
| [`Practical4_TEtranscripts.md`](Practical4_TEtranscripts.md) | Student handout. |
| `Practical4_TEtranscripts.html` | The same handout as a single HTML file. [Open it in a browser](https://clemgoub.github.io/Physalia_TE_course/Practical4/Practical4_TEtranscripts.html), or download it to read it offline. |
| [`INSTRUCTOR_practical4.md`](INSTRUCTOR_practical4.md), [`.html`](https://clemgoub.github.io/Physalia_TE_course/Practical4/INSTRUCTOR_practical4.html) | Instructor run sheet: the commands to type live, and how long each took on the course server. |
| [`VolcanoPiwi.R`](VolcanoPiwi.R) | The `R` script that draws the volcano plots. Students copy it from the course server. |
| [`figures/`](figures/) | The volcano plot of the handout, drawn with `VolcanoPiwi.R`. |

To redraw the volcano plot, from a directory that has the results in `results/`:

```bash
Rscript -e 'png("figures/volcano_piwiKD.png", width = 1400, height = 700, res = 110); source("VolcanoPiwi.R"); dev.off()'
```

To rebuild the HTML, from this directory:

```bash
../render/render.sh Practical4_TEtranscripts.md "Practical 4: TEtranscripts" Practical4_TEtranscripts.html
```
