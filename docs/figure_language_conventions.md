# Figure Language Conventions

These conventions apply to notebook-aligned recovery figures intended for dissertation use.

## Core Rule

Figure-internal text should contain only necessary scientific labels. Interpretive, explanatory, workflow, or tutorial-style wording belongs in captions or surrounding text, not inside the figure.

## Titles

Use concise scientific titles:

- `Thomas-Fermi condensate`
- `Scalar phase map`
- `PCI image`
- `DGI image`
- `Faraday dark-field image`
- `Faraday dual-port signal`
- `Camera image`
- `Noisy PCI camera frame`
- `Multishot sequence evolution`

Avoid:

- `Stage X`
- `Step X`
- `Notebook-aligned recovery`
- `reference implementation`
- explanatory subtitles such as `what every imaging mode actually integrates over`

## Panel Titles

Panel titles should identify the quantity, not explain the workflow.

Preferred examples:

- `(a) 3D density cuts`
- `(b) x-integrated column density`
- `Scalar phase map`
- `Central lineouts`
- `Faraday dual-port signal`

## Legends

Use concise scientific labels and notation:

- `$x$ cut`
- `$y$ cut`
- `$z$ cut`
- `$R_x=1.19\,\mu\mathrm{m}$`
- `$R_y=24.82\,\mu\mathrm{m}$`
- `$R_z=1.49\,\mu\mathrm{m}$`

Avoid loose or colloquial labels such as `what we see`, `actually`, `across`, or `along` unless those words are needed in a caption.

## Notation

Use Matplotlib mathtext-compatible notation:

- full column-density distributions: `$n_{\mathrm{col}}(y,z)$`, `$n_{\mathrm{col}}(x,z)$`, `$n_{\mathrm{col}}(x,y)$`
- peak column-density scalars: `$\tilde{n}_x$`, `$\tilde{n}_y$`, `$\tilde{n}_z$`
- scalar phase: `$\phi$`
- Faraday angle: `$\theta_F$`
- normalised intensity: `$I/I_0$`
- atom number: `$N_0$`
- pulse duration: `$\tau$`

Do not use peak-scalar notation such as `$\tilde{n}_x$` for a full 2D column-density map.

## Units

Use mathtext for units:

- `$\mu\mathrm{m}$`
- `$\mathrm{m}^{-2}$`
- `$\mathrm{m}^{-3}$`
- `$\mathrm{cm}^{-2}$`
- `$\mathrm{cm}^{-3}$`
- `$\mu$s`
- `$e^-$`

Avoid plain-text superscripts such as `m^-2`, `m^(-2)`, or mixed `um` / Unicode `μm` labels in figure text.

## Number Formatting

Displayed numerical values should be rounded for readability:

- radii: usually two decimal places in legends
- frame quantities: use compact labels such as `25.0k` atoms and integer nK temperatures
- avoid long decimals in titles or legends

Detailed numerical values belong in CSV, JSON, tables, or captions.

## Captions

Captions should carry interpretation:

- which notebook stage or script generated the figure
- whether a figure is an exact recovery or a model extension
- whether outputs are calibrated or representative
- scientific interpretation and limitations

Do not put those statements inside the figure unless they are a short scientific label.
