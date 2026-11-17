# Sandframe: Game Off 2026 project brief

Working title "Sandframe" (not cleared, see Open questions). This file is the single brief for the whole jam, from the physics toy to the submitted game. Claude Code reads it at the start of every session.

## How to work with this file

- The user, Chris, is a fluid-mechanics professor and a solo part-time developer. He cannot run lab experiments; he judges the physics by watching renders and PNGs. Make the physics visible and every parameter easy to change.
- Work only on the phase marked CURRENT in Status. Do not start the next phase without the user's go.
- At the end of each task: run the tests, then report what passes, every assumption you made that is not in this file, and anything here you think is wrong.
- Never change an entry in Decisions on your own; propose the change and wait.
- At the end of each session, append one dated line to Status.
- Do not add features that are not in the current phase. When in doubt, ask.

## Status

- CURRENT: Phase 0
- 2026-11-17: brief written; nothing built yet.
- 2026-11-17: Phase 0 built (repo celoy/sandframe, .venv on Python 3.11.7, pytest smoke test passes, empty Godot 4.7.2 project exports to build/web); waiting for the user to open the web export in Chrome.

## Deadline and jam

- GitHub Game Off 2026, theme SHIFT. Submissions close 1 December 2026 at 21:37 UTC (13:37 PST). Target submission: 30 November, evening (Paris time), leaving a day of buffer.
- The source must be in a public GitHub repository. Any engine is allowed.
- Rated categories: overall, gameplay, graphics, audio, innovation, theme interpretation.
- Before submitting, the user checks the jam page for any rule on AI-assisted code (open question).

## The game in one page

Pitch. You hold a moving-sand frame: layers of coloured sand in a viscous liquid between two glass plates. Shift it and flip it so that the sand settles into the landscape sketched on the glass.

Verbs (only two):
- Shift: in-plane rotation of the frame, continuous, controlled by the player, rate-limited, capped at ±60°.
- Flip: 180° rotation about the horizontal axis in the plane of the frame. Instantaneous, physics frozen during it, shown as a 0.3 s vertical squash. Inverts up and down; left and right are kept.

One level:
1. The level gives a starting deposit (layered sand on the floor) and a pencil sketch of the target outline.
2. The player flips. The deposit becomes a reservoir hanging from the ceiling, which releases in plumes. The player shifts during the fall to aim the curtains, and before the next flip to reshape the deposit (avalanches at the repose angle; past it, the whole landscape pours towards a wall).
3. Repeat. The level is judged when the sand has settled and the shift is back to 0°.

Scoring:
- Outline: the settled surface h(x) is compared with the target h*(x), allowing a horizontal slide s within ±slide_max. Score = max over s of 1 − Σ|h(x + s) − h*(x)| / Σ h*(x). The stripes inside are not judged.
- Drapes: each flip ends with a thin fine-sand drape across the frame. The sketch marks summits with the number of white drapes required.
- Stars: outline ≥ threshold; required drapes present; flips ≤ par. Thresholds are set from bot runs.
- Free mode: no target; flip, shift, watch, save a postcard PNG. This is the contemplative side.

Why levels are hard (design levers, derived from the model):
- A shift during a fall translates the whole curtain pattern; it cannot change the number of peaks. Gathering sand needs a pour past repose, then a flip.
- The reservoir empties last-in, first-out, so what lands last is decided one cycle earlier.
- One drape per flip, so the number of drapes on a summit is a visible lower bound on flips.
- Mass is conserved: a peak larger than any region of the reservoir needs gathering.
- If plumes are seeded by the former summits (to be checked in the toy), the landscape you build before a flip is the program of the next release.

Fun hypothesis (falsifiable with bots, Phase 3):
- A search bot over shift sequences beats a flip-only bot by a clear margin (target: 0.15 in outline score) on most levels.
- For a fixed action sequence, the spread of scores across random seeds is smaller than the gap between star thresholds.
- On the hard levels, a greedy bot (one step of lookahead) fails where the search succeeds.
- The set of outlines reachable with 2–3 flips is varied (no two generated levels look alike).

Kill criterion and fallback: if by 25 November the game model is not controllable (search bot within 0.05 of the flip-only bot), drop the levels and submit the contemplative sandbox: free mode, shift, flip, postcard export.

Out of scope for the jam: full in-plane rotation beyond ±60°, bubbles, more than three sand classes, level editor, menus beyond level select, settings screen, voice, 3D, save export and import, mobile touch controls.

## Physics reference model (the toy)

Python toy that shows the physics; it is the reference against which the game model is calibrated. Clarity and correctness over speed.

System: closed cell of width W and height H in the frame plane, gap b between the plates, liquid of density rho_f and viscosity mu_f, three sand classes i with diameter d_i and density rho_i, volume fractions phi_i(x, y, t), phi = Σ phi_i.

Assumptions (repeat them in the docstrings):
1. Gap-averaged Darcy flow (Hele-Shaw), permeability K = b²/12. Valid for structures much larger than b and negligible inertia.
2. The suspension is a Newtonian mixture of constant effective viscosity mu_m = mu_r · mu_f. Real viscosity diverges near packing and dense layers can have a yield stress; ignored.
3. Boussinesq: density differences only in buoyancy.
4. Relative to the mixture, grains settle along gravity with hindered Stokes velocity V_i (1 − phi)^n, V_i = (rho_i − rho_f) g d_i² / (18 mu_f), n = 4.65 (Richardson–Zaki, recalled). Particle Reynolds numbers stay below 0.1 with the defaults.
5. Isotropic diffusion D regularises the finger scale. For a sharp interface in Darcy flow with equal viscosities, n(k) ≈ K g Δρ k / (2 mu_m) − D k² (derived), maximal at lambda_star when D = lambda_star K g Δρ / (8 π mu_m). lambda_star is a guess (default 3b); D is computed from it.
6. The deposit (bed) is rigid, impermeable, packed at phi_bed; no compaction.
7. The flip is instantaneous; nothing moves relative to the frame during it.

Frame coordinates: x along the frame, y up the frame. Shift angle theta > 0 turns the frame anticlockwise on screen. Gravity in frame coordinates: g_vec = (−g sin theta, −g cos theta).

Flow. Density anomaly of the suspension, bed excluded: Δρ = Σ phi_i (rho_i − rho_f). Darcy: u = −(K/mu_m)(∇p − Δρ g_vec), ∇·u = 0. With u = ∂ψ/∂y, v = −∂ψ/∂x:

    ∇²ψ = (K/mu_m) (g_x ∂Δρ/∂y − g_y ∂Δρ/∂x),   ψ = 0 on the walls and in the bed.

At theta = 0 this is ∇²ψ = (K g/mu_m) ∂Δρ/∂x. Sign check: a heavy region on the right sinks on the right, with clockwise circulation. Solver: 5-point Laplacian on fluid cells only (scipy.sparse), `splu` factorisation reused until the bed mask changes; a DST solver on the full rectangle is fine before the first deposit.

Transport, each class: ∂phi_i/∂t + ∇·(phi_i (u + w_i)) = ∇·(D ∇phi_i), with w_i = V_i (1 − phi)^n g_vec/g. Finite volume, conservative, upwind with minmod limiter, explicit, CFL ≤ 0.5 and diffusion number ≤ 0.25. Zero flux through walls. Keep phi_i ≥ 0 and phi ≤ phi_bed (limit inflow into nearly full cells; document how).

Deposition: settling flux through any face between a suspension cell and a bed cell leaves the suspension and is added to the bed column it enters, as thickness flux / phi_bed. A suspension cell touching the bed with phi ≥ 0.95 phi_bed is converted entirely into a bed layer. No advective flux enters the bed (ψ = 0), so dense plumes spread along the deposit as gravity currents; this is wanted.

Bed: one stack of layers per column (composition fractions, thickness) and a surface height h(x), measured along y. Cells whose centre lies below h(x) are bed cells. Merge adjacent layers with nearly identical compositions. Avalanches: with surface angle alpha = atan(dh/dx), the surface is stable if |alpha + theta| ≤ theta_r. While a pair of neighbouring columns violates it, move volume from the top layers of the higher column to the lower one, top layers first. Past theta_r, the whole landscape pours towards the low wall; the column representation holds within ±60°.

Flip: y → H − y for all suspension fields. Each bed column becomes suspension hanging from the ceiling at phi_bed, layers in reversed order (oldest against the ceiling). The bed becomes empty. Add a seeded random perturbation of relative amplitude interface_noise in the interface cells. A flip may happen at any time. Shift angle is unchanged by a flip.

Toy parameters (params/default.toml; mark each as derived or guess in comments):

| Name | Default | Status |
|---|---|---|
| W, H | 0.20 m, 0.10 m | guess |
| nx, ny | 200, 100 (dx = 1 mm) | guess |
| b | 3 mm | guess |
| rho_f, mu_f | 1150 kg/m³, 0.02 Pa·s (glycerol–water) | guess |
| coarse | d = 300 µm, rho = 2650, dark brown | guess |
| medium | d = 150 µm, rho = 2650, coral | guess |
| fine | d = 50 µm, rho = 2500, white | guess |
| mu_r | 30 | guess (literature fits up to about 100 in water) |
| phi_bed | 0.55 | guess |
| n_rz | 4.65 | recalled |
| theta_r | 30° | guess (underwater repose) |
| lambda_star | 3 b | guess; D derived |
| interface_noise | 0.05 | guess |
| omega_max | 60°/s shift rate | guess |
| seed | 1 | |

Print at start-up: K, D, plume velocity scale U = K g Δρ / mu_m (about 1 cm/s), each V_i and fall time over H (about 27 s, 110 s, 16 min), particle Reynolds numbers, time step, linear growth rate at lambda_star.

Toy scenarios: flat_bands (three horizontal bands, coarse at the bottom), two_hills (layered hills of different heights; tests whether former summits seed the first plumes), uniform_suspension (sorting and final drape), random_hills (seeded).

Toy viewer (pygame-ce): frame at about 5 px/mm drawn rotated by theta, side panel. Suspension colour mixes class colours over the liquid colour weighted by phi_i; the bed is drawn from the column layers so stratigraphy stays sharp. Keys: F flip (optional squash), Q/E or ←/→ shift (rate-limited), 0 return shift to 0°, Space pause, +/− speed, 1–4 scenarios, R reset, N new seed, V velocity arrows, C total-phi heat map, I interface contour, P save PNG, L reload parameters. Panel: time, speed factor, shift angle, flips, relative mass error per class, suspended fraction per class, dominant interface wavelength.

Toy headless runner: `python -m sandframe_toy.cli --scenario two_hills --actions "flip;wait 30;shift 40 10;shift 0 5;wait 600" --snap-every 5 --out runs/two_hills_001` writes PNG snapshots and diagnostics.csv (t, mass per class, suspended fraction per class, dominant wavelength, plume front position and speed, max |u|, shift angle).

Toy tests (pytest):
1. Poisson solver against an analytic sine mode, with and without bed mask.
2. Sign test at theta = 0 and at theta = 45°.
3. Mass of each class conserved to 1e-9 relative over 1000 steps including deposition, avalanches, shifts and a flip.
4. phi_i ≥ 0 and phi ≤ phi_bed always.
5. Two flips with no time between them return the initial state.
6. Avalanches conserve volume; at rest all surfaces satisfy |alpha + theta| ≤ theta_r.
7. Linear growth of a small sinusoidal interface perturbation within 30 % of n(k) for resolved wavenumbers.
8. Settling only: a uniform single-class suspension descends at the hindered velocity.

## Game model (decided at Gate A)

The toy is too heavy for GDScript in a single-threaded web export. The game runs one of two models, chosen at Gate A. Whichever is chosen, it is first written in Python (`model/`, numpy only) as the twin used by the level generator and the bots, then ported to GDScript and checked against golden files produced by the twin.

Route A, coarse Darcy port. The toy equations on a grid of about 100 × 50, Poisson solved by conjugate gradient or red-black SOR on the masked domain, in GDScript. Chosen only if the benchmark meets the budget.

Route B, reduced model (default). Columns j = 1..N (N about 120) along the frame.
- Bed: stack of layers per column, as in the toy.
- Reservoir after a flip: the bed stacks hanging from the ceiling, reversed.
- Plume sites: local downward bulges of the reservoir's lower interface with prominence > p_min, plus background sites about lambda apart with seeded jitter where no bulge lies within lambda; sites closer than lambda_merge merge.
- Release: each site drains its catchment (columns between midpoints to neighbouring sites) at rate Q, taking the lowest reservoir layers first.
- Split: of the released class i, a fraction f_i falls in the plume as parcels at speed U_p along gravity; the rest joins a dilute cloud (amount per column per class) that settles at V_i along gravity.
- Landing: parcels and cloud drift sideways by (fall distance) · tan theta; parcels deposit with a spreading kernel of width w_s (gravity current), then avalanches act as in the toy.
- Recirculation is not modelled; its effects are absorbed in lambda, f_i and w_s.
Calibrated parameters: lambda (vs b), lambda_merge, p_min, Q, U_p, f_i per class, w_s, drape thickness per flip. Each comes from a named toy experiment in tools/ with a CSV and a fitting script.

Gate A (planned 22 November). Inputs: the user's approval of the toy physics; the calibration fits; a GDScript benchmark of one coarse Darcy step in a web export on the user's laptop. Budget (guess): ≤ 6 ms per sim step, with 2 steps per rendered frame. Route A if the budget is met and the user wants it; otherwise Route B. Also check at Gate A whether plumes start under former summits; if not, revise the design levers in "The game in one page" before Phase 3.

## Levels, scoring and generator

Level file (game/levels/Lxx.json):

```
{
  "id": "L03",
  "seed": 3,
  "sands": ["coarse", "medium", "fine"],
  "initial_bed": [ [ {"comp": [0.7, 0.3, 0.0], "t": 0.004}, ... ], ... ],
  "target": { "h": [ ... N values in m ... ], "drapes": [ {"x": 0.06, "count": 2} ] },
  "par_flips": 2,
  "slide_max": 0.02,
  "thresholds": { "outline": 0.80 },
  "solution": [ "flip", "shift 40 3.0", "shift 0 2.0", "wait 60", "flip", "wait 120" ],
  "meta": { "greedy_fails": true, "min_flips_found": 2 }
}
```

Drape detection: in a column's stack, a layer with fine fraction > 0.6 and thickness > t_min counts as a drape.

Generator (tools/generate_levels.py, on the Python twin):
- Action alphabet: shift to one of {−60, −30, 0, 30, 60}° and hold {short, long}; wait; flip. Shifts are rate-limited as in the game.
- Candidate targets: settled outlines of random action sequences, so every level is solvable.
- Par: smallest number of flips found by beam search (width about 200) that reaches the outline threshold. This is an upper bound; keep a target as par k only if the best score found with fewer flips is below threshold minus delta (delta = 0.05, guess).
- Hard levels: those where the greedy bot fails.
- Report per level: par, greedy result, score spread over 10 seeds for the stored solution.
- Produce about 15 candidates; the user picks 8–10 and orders them.

Bot runs for the fun hypothesis: flip-only bot, greedy bot, beam search, random bot; write a CSV and a short summary table.

## Godot game

- Godot 4.3 or later (confirm the user's version), GDScript only, Compatibility renderer, web export with thread support off. Flag any Godot 3 idiom (yield, export without @, onready without @, PoolByteArray, KinematicBody...).
- Simulation in plain classes (RefCounted), no Nodes, testable headless: Params, SandState, GameModel (route A or B), Scoring, LevelLoader. Scenes only display state and pass input.
- Scenes: Title, LevelSelect, Play, Free. Minimal text everywhere.
- Frame view: the sim state rasterised into an Image via PackedByteArray and an ImageTexture updated each frame; the frame node rotates for the shift; the flip is a scale.y tween through zero, physics paused.
- Target sketch: dashed pencil Line2D in frame coordinates; it rotates with the shift but is not flipped (open question: overlay vs glass).
- Input: drag the frame rim, mouse wheel or A/D and ←/→ to shift (rate-limited, omega_max); F flip; Space pause; Tab speed ×1/×4/×16; Z undo last flip (state snapshot before each flip); R restart; Enter judge (only when settled and shift is 0°, otherwise the frame returns to 0° first).
- HUD: flips used against par as dots, stars, speed. Free mode: P saves a postcard PNG.
- Audio: sand hiss whose volume follows the sand flux (plumes plus deposition). Try AudioStreamGenerator; if it misbehaves in the web export, loop a CC0 sample and modulate its volume. Credit every CC0 asset in the repo.
- Progress: stars per level in user://progress.json.
- Performance budget: sim ≤ 8 ms per frame on average in Chrome on the user's laptop (guess); test Chrome and Firefox.
- Tests: `godot --headless --script res://tests/run_tests.gd` runs unit tests and compares the GDScript model with golden JSON files from the Python twin (same initial state and actions, tolerance documented).

## Repository layout

```
sandframe/
  CLAUDE.md
  README.md            # written in Phase 6
  LICENSE              # chosen by the user
  toy/                 # Python reference physics (package sandframe_toy)
    params/default.toml
    src/sandframe_toy/ params.py state.py flow.py transport.py bed.py flip.py shift.py
                       scenarios.py diagnostics.py sim.py cli.py
    app/view.py
    tests/
  model/               # Python twin of the game model (package sandmodel), numpy only
    src/sandmodel/
    tests/
  tools/               # calibrate_*.py, generate_levels.py, bots.py, golden.py
  game/                # Godot project
    project.godot
    sim/  scenes/  ui/  levels/  audio/  tests/  golden/
  docs/                # screenshots and GIFs for the itch page
  runs/                # gitignored
```

## Conventions

- Python 3.11+, one virtual environment at the repo root (.venv). Dependencies: numpy, scipy (toy only), pygame-ce, pytest. Nothing else without asking.
- SI units inside every model; units in names and docstrings.
- Parameters in TOML (Python) and JSON or .tres (Godot), never hard-coded.
- One seed per run; every run reproducible.
- Vectorise with numpy; Python loops over columns are fine for avalanches.
- Content (levels, texts) in JSON.
- Commit at the end of each task with a message naming the phase and task.

## Plan

Dates and hours are guesses. If a phase runs more than a day late, apply the cut list.

Phase 0, setup (Tue 17 Nov, about 1 h)
- Public GitHub repo, layout above, .gitignore, .venv, pytest running on an empty test, empty Godot project that exports to web and runs in Chrome.
- Acceptance: the user opens the empty web export locally.

Phase 1, toy flip (Wed 18 – Fri 20 Nov, about 6 h)
- Step 1: params, state, flow, transport; toy tests 1, 2 (theta = 0), 3, 4, 7, 8.
- Step 2: bed, deposition, avalanches, flip, scenarios; tests 3, 5, 6.
- Step 3: viewer. Step 4: headless runner and diagnostics.
- Acceptance: tests pass; the user reviews renders of the four scenarios and his observations are recorded in Status.

Phase 2, toy shift, calibration, benchmark (Sat 21 – Sun 22 Nov, about 5 h)
- Shift in the toy (gravity rotation, tilted avalanches, rate limit); test 2 at 45°, tests 3 and 6 with shifts.
- Calibration experiments with CSV outputs and fits: finger spacing vs b; merge distance; plume seeding at bulges (two_hills); plume speed U_p; release rate Q; class fractions f_i in plumes vs dilute cloud; spreading width w_s; drape thickness per flip.
- GDScript benchmark of one coarse Darcy step in a web export.
- Gate A on Sunday 22: the user chooses Route A or B; record the decision.

Phase 3, game model twin, bots, levels (Mon 23 – Wed 25 Nov, about 5 h)
- model/ implementing the chosen route, with tests (mass conservation, flip twice = identity, avalanche rest state, drape count = flips on a simple case).
- Comparison plots twin vs toy on the four scenarios.
- Bots and generator; 15 candidate levels; fun-hypothesis summary table.
- Golden files for the Godot tests.
- Kill check on Wed 25 (see Kill criterion).

Phase 4, Godot sim and view (Wed 25 – Fri 27 Nov, about 6 h)
- GDScript port of the model; golden tests pass; frame view, shift input, flip squash, pause and speed.
- Acceptance: free mode playable in a web export within the performance budget.

Phase 5, game loop (Fri 27 – Sat 28 Nov, about 4 h)
- Level loader, target sketch, judging, stars, par dots, undo, restart, level select, progress save, audio.
- Acceptance: all chosen levels finishable in the web export; the solution of each level, replayed, earns its stars.

Phase 6, freeze, export, submit (Sun 29 – Mon 30 Nov, about 3 h; Tue 1 Dec buffer)
- Feature freeze on Sunday 29. Bug fixes only after that.
- Two or three outside players try it without help; note where they get stuck.
- README (how to play, physics in one paragraph, credits), itch page text and screenshots, web build uploaded, submission made on Monday 30.

Cut list, in this order: undo; AudioStreamGenerator (use a looped sample); level generator (5 hand-made levels with par found by the beam search only); drape star; route A (fall back to B); levels altogether (fallback: contemplative sandbox).

## Decisions

- Godot Wild Jam #99 dropped (no AI-written code allowed); entry for GitHub Game Off 2026, theme SHIFT.
- Kind of game: indirect control in real time with pause; the player shapes a landscape; a natural process with a visible plan moves across it; judged when it comes to rest; contemplative, soft score.
- Moving-sand frame in liquid chosen; sand poured in air and rain-erosion versions dropped.
- Target is the outline only, matched up to a horizontal slide; stripes are not judged.
- Shift: in-plane rotation, continuous, player-controlled, rate-limited, capped at ±60°.
- Flip: 180° about the horizontal in-plane axis, instantaneous, physics frozen, vertical squash animation; up/down inverted, left/right kept.
- Par is the number of flips; drapes (one per flip) give a visible lower bound; the solver gives the upper bound.
- Physics hypotheses to check in the toy: plumes seeded by the downward bulges of the inverted deposit; two-phase fall (collective plumes, then dilute settling with fines last as a frame-wide drape).
- The toy uses a gap-averaged Darcy mixture model as reference physics.
- The game model is written first as a Python twin, then ported to GDScript and checked against golden files.

## Open questions (as tasks)

- Target sketch: overlay attached to the frame (default) or drawn on the glass and flipped with it. Decide after first playtest of Phase 5.
- Title: search Steam, itch.io and the web for clashes before adopting "Sandframe" or any other title.
- Licence for the public repo.
- Confirm the jam's rule on AI-assisted code before submitting.
- Godot version installed.
- Route A or B (Gate A).
- Star thresholds and delta (from bot runs).
- Audio approach in the web export.

## Sources

- A. Lange, M. Schröter, M. A. Scherer, A. Engel, I. Rehberg, "Fingering instability in a water-sand mixture", arXiv:patt-sol/9706004, https://arxiv.org/abs/patt-sol/9706004 (opened). Flip of a sand layer in a Hele-Shaw-like cell; growth rates rise then plateau with k; effective viscosity fit about 100 times water; weak dependence on grain diameter; dense packings can have a yield stress.
- C. Völtz, W. Pesch, I. Rehberg, "Rayleigh-Taylor instability in a sedimenting suspension", Phys. Rev. E 65, 011404 (2001), https://epub.uni-bayreuth.de/4031/ (opened). Same geometry with glycerin.
- H. A. Makse, S. Havlin, P. R. King, H. E. Stanley, "Spontaneous stratification in granular mixtures", Nature 386, 379 (1997), https://doi.org/10.1038/386379a0 (seen in search results). Stratification when pouring in air; cosmetic here.
- Game Off 2026 page, https://itch.io/jam/game-off-2026 (seen in search results): submission window.
- Commercial moving-sand frames: flow times of 15 min to 2 h (seen in search results).
- Recalled, not verified: Richardson–Zaki hindered settling (1954); Engelstein and Shalev, Building Blocks of Tabletop Game Design (2019), input vs output randomness; Koster, A Theory of Fun (2004).
