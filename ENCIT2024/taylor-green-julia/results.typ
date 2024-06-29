#let ivory = rgb("#e6deca")
#set page(paper: "a4", fill: ivory)
#set text(font: "Crimson Pro", size: 11pt)
#set par(justify: true, leading: 0.65em)

#show raw: txt => text(fill: rgb("#000048"), weight: 550, txt)

= Computing Environment Description
#v(6pt)

Host `freebird` is a
`Gigabyte` computer `desktop`; with
`Intel Core i3-9100F` CPU model:
`64` bits,
`core 9` gen,
`1` CPU with
`4` cores with
`256KiB` L1 cache (`4x32 KiB` for data and `4x32 KiB` for instructions),
`1024 KiB` L2 cache, and
`6 MiB` L3 cache, running at
`3299 MHz` (avg.); with
`16 GiB` total RAM (
`15.57 GiB` available, and
`4.64 GiB` used by some `234` concurrent processes) on a
`Linux` operating system with a
`6.6.32-1-MANJARO` kernel,
`x86_64` architecture, compiled with
`gcc` version
`14.1.1` with
`60% (default)` of swappiness, with
`16 GiB` of hard disk swap of which
`0.0 MiB (0.0%)` were in use during benchmark tests.

Host `stilo` is a
`Positivo Informatica SA` computer `laptop`; with
`Intel Celeron N3010` CPU model:
`64` bits,
(unknown) gen,
`1` CPU with
`2` cores with
`112KiB` L1 cache (`2x24 KiB` for data and `2x32 KiB` for instructions),
`2 MiB` L2 cache (`2x1024 KiB`), and
no L3 cache, running at
`539 MHz` (avg.); with
`4 GiB` total RAM (
`3.69 GiB` available (est.), and
`2.65 GiB` used by some `213` concurrent processes) on a
`Linux` operating system with a
`5.15.130-1-MANJARO` kernel,
`x86_64` architecture, compiled with
`gcc` version
`13.2.1` with
`60% (default)` of swappiness, with
`4 GiB` of hard disk swap of which
`518.9 MiB (12.7%)` were in use during benchmark tests.

= Benchmarks on 2D Taylor-Green Vortex Decay with `D2Q9` LBM

== C (ISO C99) Codes

The `ISO-C99` code was compiled with `GCC`, version `14.1.1` of `2024-05-22` (Copyright #sym.copyright 2024 Free Software Foundation,  Inc.)  in
each one of the computing environments described above. ...

== Julia Codes

The `julia` language tests were  performed  with  the  following  characteristics:  version  `1.10.2`  (commit  `bd47eca2c8a`  official  build),
`libopenlibm`, and `libLLVM-15.0.7`, reporting, respectively on hosts `stilo`, and `freebird`, CPU frequencies of `1.04GHz`, and of `3.60GHz`.

#figure(
  table(align: center + horizon,
    columns: 8,
    inset: 4pt,
    table.header[*Domain size*][*Time Steps*][*Precision*][*Memory*][*Code*][*Host*][*Speed (Mlups)*][*Relative Speed*],
    table.cell(rowspan: 4)[32 #sym.times 32],
    table.cell(rowspan: 4)[204'800],
    table.cell(rowspan: 4)[`Float64`],
    table.cell(rowspan: 4)[168.47 KiB],
    [`OP1`], [`stilo`]                         , [ 3.61], [1.000],
    [`ORI`], table.cell(rowspan: 3)[`freebird`], [21.69], [0.916],
    [`OP1`],                                     [23.68], [1.000],
    [`OP2`],                                     [20.72], [0.875],
  ),
  caption: [Benchmark results for tested julia codes],
) <julia-times-1>

= Implementation Lessons from the Benchmarks

== Julia

The `OP1` code optimization set include:

- Explicitly declaring a global `const` precision type, named `𝕋`, and using it in all type  declarations  and  whenever  initializing  floating
  point quantities, as in `const 𝕋 = ℙ` --- with `ℙ` being previously defined in the REPL session, as a floating point precision type, as in
  `ℙ = Float64`.
- Replacing hardcoded constants by a `const`, typed variable consistent with its usage, as in `const chunk = UInt(32)`.
- Saving results to a local variable instead of calculating it repeated (=`ndir`) times inside double loops, as in `𝘂𝘂 = 𝚞 * 𝚞 + 𝚟 * 𝚟`  outside
  the `for` loop that uses it. This might be the greatest source of increased speed.
- Using `tuple`-like multiple assignments whenever convenient, as in `ϱ, 𝚞, 𝚟 = ρ[i], 𝑢[i], 𝑣[i]`.
- Using cascading initialization `=` assignments whenever convenient, as in `ϱ = 𝚞 = 𝚟 = zero(𝕋)`.

It's implementation resulted in a gain of relative speed from $0.916$ to $1$, for the $32 #sym.times 32$, `Float64` case, as shown on
@julia-times-1.

Further attempts at optimization focused mainly on declaration placements where gathered in a set named `OP2`. As the benchmarks  reveal,  `OP2`
resulted in such drastic performance _loss_ that outweighted the gains earned with `OP1`. Therefore, the strategy implemented on `OP2` is to  be
avoided:

// !j 144 -i2 -H-2
- Moving some function's unique `local` constants into the global scope, even if they are declared as  `const`s,  and  being  explicitly  typed,
  access to a function's locals is still faster, even if they have to be re-computed time and again when their defining function  is  repeatedly
  called.


