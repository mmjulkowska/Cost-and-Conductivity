# Figures

## GRANAR anatomies
`plot_granar.py` draws every GRANAR anatomy, with cells coloured by tissue and with exodermis and endodermis outlined wherever MECHA assigns a barrier to that segment. It writes:
- `outputs/figures/granar/<id>.png`: one image per anatomy (208).
- `outputs/figures/granar_gallery_<segment>.pdf/.png`: every replicate of a segment on one page, all at the same µm scale so root sizes can be compared.

```
python figures/plot_granar.py              # all; or: python figures/plot_granar.py MP1_cell1
```

## Figure icons

`make_icons.py` draws one icon per segment: the left half is one of Emily Miller's photographs, the right half is the GRANAR anatomy built from that same section, with barrier walls coloured where MECHA assigns them.
- Each section was chosen for clarity and for a root diameter near its segment's median.
- The replicate shown is the one whose kr is closest to the segment median.

`icon_sections.json` lists, for each section:
- the photo's path inside the Box "Cross Sections" folder;
- the root centre and radius in pixels.

The script needs the photos in `data/Cross Sections/`, plus the outputs of steps 1 and 2. It writes to `outputs/figures/icons/`.

`fig_icons.py` builds `outputs/figures/Fig_conductivity_barriers_v5_icons*`. Circle size is proportional to the photographed root diameter.

```
python figures/make_icons.py && python figures/fig_icons.py
```
