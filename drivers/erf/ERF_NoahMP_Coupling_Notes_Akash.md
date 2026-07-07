# ERF ↔ Noah-MP coupling (drivers/erf) — what we changed and why

This document describes the changes on this branch to the ERF driver of Noah-MP
(`drivers/erf/`), how they relate to the upstream effort in **erf-model/ERF PR #3405**
("Expose IO and precipitation variables from Noah-MP, enhance ERF–Noah-MP coupling for
precipitation"), and the physics reasoning behind each change.

The ERF C++ side that pairs with this driver lives in the ERF repository
(`Source/LandSurfaceModel/Noah-MP/ERF_NOAHMP.{H,cpp}`); the two must be read together
because the precipitation/flux contract spans the C↔Fortran boundary.

---

## 0. TL;DR — relation to PR #3405

PR #3405 and this branch solve the **same problem** (feed precipitation into Noah-MP and
expose its IO so the LSM actually does something across the ERF boundary). We agree with
#3405 on the most important structural point — **surface fluxes are kinematic (÷ρ)** — but
this branch additionally fixes three things #3405 does not:

| Topic | PR #3405 (erf-model, open) | This branch | Consequence of the difference |
|---|---|---|---|
| Momentum/moisture flux ÷ρ | ✅ `HFX/(ρCp)`, `LH/(ρL)`, `TAU/ρ` | ✅ same | (consistent) |
| **Exner factor on heat** | ❌ feeds `HFX/(ρCp)` = kinematic **T**-flux | ✅ `HFX/(ρCp)·(p0/p)^(Rd/Cp)` = kinematic **θ**-flux | #3405 feeds a temperature flux into ERF's *potential*-temperature energy variable → biased surface heating |
| **opt_snf=4 precip partition** | ❌ lumps: `MP_SNOW=MP_GRAUP=0` in the driver | ✅ real microphysics `MP_RAINNC/MP_SNOW/MP_GRAUP` from the caller | under `opt_snf=4`, #3405 computes `FrozenPrecipFrac=0` → **all precip falls as rain, never snow** |
| **RAINBL/SR/MP_* kind** | fields bound to C `double*` (kind not audited here) | forced `real(kind=C_DOUBLE)` explicitly | a real4 field bound as `double*` reads garbage (±1e10–1e25) → water-balance abort (we hit this) |
| Restart completeness | not in scope of #3405 | full prognostic + scheme-gated state checkpoint I/O | ERF restarts do not cold-start the LSM |

Everything below is the detail behind that table.

---

## 1. Precipitation into the LSM

### What was there before
The ERF driver fed Noah-MP **no precipitation**: `RAINBL` was never set, so the land model
ran precip-free and soil could only ever dry out over an event. `MP_RAINNC/MP_SNOW/MP_GRAUP`
were derived inside `NoahmpDriverMainMod.F90` by lumping:

```fortran
RAINNCV    = RAINBL          ! (RAINBL itself was ~0)
SNOWNCV    = SNOWBL          ! SNOWBL = 0
GRAUPELNCV = 0.0
...
MP_RAINNC  = RAINNCV
MP_SNOW    = SNOWNCV         ! -> 0
MP_GRAUP   = GRAUPELNCV      ! -> 0
```

This is also what PR #3405's driver still does (`MP_SNOW = SNOWNCV`, `MP_GRAUP =
GRAUPELNCV`, both zero).

### What WRF does (the reference)
WRF's Noah-MP *caller* (`drivers/wrf/NoahmpWRFmainMod.F90:553-568`) sets
`MP_RAINNC/MP_SNOW/MP_GRAUP/MP_RAINC/MP_SHCV/MP_HAIL` **directly from the microphysics**,
where snow and graupel are *subsets* of the total non-convective `MP_RAINNC`. WRF's
driver-main does **not** re-derive them.

### What we did
We made the ERF path do what WRF's caller does:

- The ERF C++ driver now computes, per land-call interval, the accumulated-precip deltas in
  **[mm]** from the microphysics accumulators and feeds all of them across the boundary:
  - `RAINBL   = Δrain_accum`  (total surface precip, mm; clamped to a physical ceiling)
  - `SR       = Δfrozen / Δtotal`  (frozen fraction, used by the "other precip" path)
  - `MP_RAINNC = Δrain_accum`  (total non-convective, mm)
  - `MP_SNOW   = Δsnow_accum`  (snow+ice subset, mm)
  - `MP_GRAUP  = Δgraup_accum` (graupel subset, mm)

  ERF's Morrison microphysics already carries exactly these with the same subset
  convention: `rain_accum` (total), `snow_accum` (ice+snow), `graup_accum` (graupel).
- In this driver (`NoahmpDriverMainMod.F90`) we **removed the lumping remap** so the
  caller-supplied `MP_*` are *not* overwritten. Only the channels ERF lacks
  (`MP_RAINC`, `MP_SHCV`, `MP_HAIL` — convective / shallow-convective / hail) are zeroed,
  which is identical to WRF when those OPTIONAL arguments are absent. `RAINBL` (total) is
  still set and still feeds the `opt_snf=1` path and `ACSNOW`.

### Why this matters — `opt_snf=4`
The shared Noah-MP core (`src/AtmosForcingMod.F90`, `OptRainSnowPartition==4`) computes the
frozen fraction as
```
FrozenPrecipFrac = PrecipFrozenTot / PrecipNonConvRefHeight,
PrecipFrozenTot  = MP_SNOW + MP_GRAUP + MP_HAIL   (each /DTBL)
```
If `MP_SNOW = MP_GRAUP = 0` (the lumped path), then `FrozenPrecipFrac = 0` and **every
precipitation event is treated as rain, regardless of temperature**. Feeding the real
microphysics breakdown fixes this; `opt_snf=1` (temperature-based partition) ignores `MP_*`
and is unaffected, so the change is safe for existing `opt_snf=1` runs.

### SNOWBL is a dead field (in WRF too)
`SNOWBL` is *not* the mechanism for snow. In WRF's own driver it is only initialized to
`undefined_real` (`drivers/wrf/NoahmpIOVarInitMod.F90:595`) and **never assigned or read**
into the core. We zero it defensively (an uninitialized `SNOWBL` propagated into
`MP_SNOW/PrecipSnow ~1e20` and tripped a balance abort — see §3). Snow reaches the core via
`MP_SNOW`, not `SNOWBL`.

**Proof (not just "it's initialized" — an exhaustive occurrence list).** A single init
link would prove nothing: *every* field is declared/allocated/initialized the same way. What
makes `SNOWBL` dead is that in WRF's driver those are the *only* places it appears — it is
never assigned a real value and never read. Grepping the entire Noah-MP tree at WRF's pinned
commit (`grep -rniE SNOWBL`), the WRF driver has exactly **three** hits, all non-functional:

| Location | What it does |
|---|---|
| `drivers/wrf/NoahmpIOVarType.F90:101` | type declaration |
| `drivers/wrf/NoahmpIOVarInitMod.F90:56` | allocation |
| `drivers/wrf/NoahmpIOVarInitMod.F90:595` | set to `undefined_real` |

There is **no fourth occurrence** anywhere under `drivers/wrf/` — no assignment, no read.
Two corroborating facts from the same grep:

1. **The shared core `src/` (all the physics) contains ZERO occurrences of `SNOWBL`** — the
   core literally cannot see it.
2. **The module that actually feeds forcing into the core,
   `drivers/wrf/ForcingVarInTransferMod.F90`, does not mention `SNOWBL` at all.** Snow enters
   the core there via `MP_SNOW`/`MP_GRAUP`/`MP_RAINNC` (÷`DTBL`), e.g.
   `PrecipSnowRefHeight = MP_SNOW(I,J)/DTBL` — *that* is the snow path, not `SNOWBL`.

So `SNOWBL` in WRF is allocated-and-initialized dead weight; snow reaches the core through
`MP_SNOW`. (By contrast, the pre-fix ERF driver at this same commit *did* read it —
`drivers/erf/NoahmpDriverMainMod.F90:58,63`: `SNOWBL*DTBL → SNOWNCV = SNOWBL` — which is
exactly the lumping remap we removed, and why an uninitialized `SNOWBL` produced the ~1e20
balance abort.)

Sources (`SNOWBL` declaration/alloc/init lines), in **both** upstream repos:

- **Canonical Noah-MP** (`NCAR/noahmp`, standalone dev repo):
  [`drivers/wrf/NoahmpIOVarInitMod.F90#L595`](https://github.com/NCAR/noahmp/blob/badab7b4b51710037fc87f3dbf329b6be59b1b5a/drivers/wrf/NoahmpIOVarInitMod.F90#L595),
  alloc [`#L56`](https://github.com/NCAR/noahmp/blob/badab7b4b51710037fc87f3dbf329b6be59b1b5a/drivers/wrf/NoahmpIOVarInitMod.F90#L56),
  and the real snow path
  [`ForcingVarInTransferMod.F90#L53-L56`](https://github.com/NCAR/noahmp/blob/badab7b4b51710037fc87f3dbf329b6be59b1b5a/drivers/wrf/ForcingVarInTransferMod.F90#L53-L56).
- **WRF** (`wrf-model/WRF`): vendors Noah-MP as the submodule `phys/noahmp → NCAR/noahmp`
  (WRF [`.gitmodules`](https://github.com/wrf-model/WRF/blob/master/.gitmodules)) pinned at
  `5da0b241`, so the shipped file is
  [`phys/noahmp/drivers/wrf/NoahmpIOVarInitMod.F90#L595`](https://github.com/NCAR/noahmp/blob/5da0b241e48ecfd9a2a1bd667ed554765856d589/drivers/wrf/NoahmpIOVarInitMod.F90#L595)
  (alloc [`#L56`](https://github.com/NCAR/noahmp/blob/5da0b241e48ecfd9a2a1bd667ed554765856d589/drivers/wrf/NoahmpIOVarInitMod.F90#L56);
  snow path
  [`ForcingVarInTransferMod.F90#L53-L56`](https://github.com/NCAR/noahmp/blob/5da0b241e48ecfd9a2a1bd667ed554765856d589/drivers/wrf/ForcingVarInTransferMod.F90#L53-L56)).

---

## 2. The RAINBL / SR / MP_* kind mismatch (the bug that masked everything)

`RAINBL`, `SR`, and `MP_RAINNC/MP_SNOW/MP_GRAUP` cross the C↔Fortran boundary as C `double*`
(`NoahmpIO.H`), bound via `C_LOC` in `NoahmpIO_fi.F90`. In this build `kind_noahmp` is
**single precision** (DOUBLE_PREC undefined), so a field declared `real(kind=kind_noahmp)`
is 4 bytes. Binding a 4-byte Fortran array to a `double*` reinterprets the bytes → Noah-MP
read **garbage RAINBL (±1e10…1e25)**, which then failed the water-balance check with values
like 1e18–1e25.

Fix: every field exposed to C++ is declared `real(kind=C_DOUBLE)` in `NoahmpIOVarType.F90`
(`RAINBL`, `SR`, `MP_RAINNC`, `MP_SNOW`, `MP_GRAUP`). This is the single most important
correctness fix on the branch — without it, no precip value is trustworthy. Any future field
added to the C boundary **must** be `C_DOUBLE`.

---

## 3. Surface-flux units returned to ERF (kinematic θ-flux + Exner)

This concerns the C++ side (`ERF_NOAHMP.cpp`) but is documented here because it is part of
the same coupling contract and is where we differ from #3405.

Noah-MP returns `HFX`, `LH` [W m⁻²] and `TAU_EW/NS` [N m⁻²] (ρ included). ERF's dynamical
core advances `RhoTheta` (ρ·θ) and consumes surface fluxes in the **kinematic MOST**
convention — potential temperature, **no ρ** — verified by direct code read of ERF:

- `ERF_MOSTStress.H` documents `surf_temp_flux` as `<θ'w'>` [K·m/s] and sets
  `t_star = −surf_temp_flux/u_star` (no ρ).
- The PBL schemes (`ERF_ComputeDiffusivityMYNN25/MYNNEDMF/MRF`, `ERF_PBLHeightYSU_MRF.H`)
  form the Obukhov length `L = −θ u*³/(κ g <θ'w'>)`, which is only dimensionally correct if
  the heat flux is the **kinematic θ-flux**.
- The LSM→surface-layer bridge (`ERF_SurfaceLayer.cpp`) sets `u_star = sqrt(tau13,tau23)`,
  so `tau` must be `u*²` [m²/s²] = **TAU/ρ**.

Correct conversion (this branch):
```
t_flux = HFX/(ρ·Cp_d) · (p0/p)^(Rd/Cp_d)   // kinematic θ-flux [K·m/s]
q_flux = LH /(ρ·L_v)                        // [kg/kg·m/s]
tau13  = TAU_EW/ρ ,  tau23 = TAU_NS/ρ       // [m²/s²]
```

**vs PR #3405:** #3405 has the ÷ρ right but writes `t_flux = HFX/(ρCp)` with **no Exner** —
i.e. it feeds a kinematic *temperature* flux where ERF's energy variable is *potential*
temperature. The missing Exner factor `(p0/p)^(Rd/Cp)` biases the surface heating (≈0.8%
near the surface, larger at lower pressure). We add Exner on the heat flux only.

> Correction history (kept for auditability): an intermediate version of this work
> mistakenly *removed* the ÷ρ (believing ERF wanted ρ-weighted fluxes, reasoning only from
> the diffusion operator's interior `−ρα∂θ/∂z`). Re-reading the PBL/MOST consumers showed
> ERF is kinematic, so ÷ρ stays; the only genuine bug in the original upstream flux code was
> the missing Exner. Net change vs original upstream = **add Exner on heat only**.

This flux convention is **not** "because WRF does it" — WRF's PBL uses its own internal
convention. The reference here is ERF's own surface-layer/PBL discretization.

---

## 4. Restart / checkpoint completeness

ERF checkpoints previously persisted only a subset of Noah-MP state, so a restart
cold-started the rest (different soil/canopy/snow trajectory than a continuous run). We added
full prognostic-state checkpoint I/O (`NoahmpWriteRestartMod.F90`, `NoahmpReadRestartMod.F90`)
plus the balance accumulators (`ACC_*XY`) and the integer `ITIMESTEP`, all guarded by
`allocated()` so the restart works for **any** namelist option set (scheme-gated fields are
only written/read when their scheme is active). This gives bitwise-consistent restart of the
land model. (Related to ERF #3255.)

Also fixed here: `NoahmpWriteLandMod.F90` land-output filename width overflow (the timestamp
string was `I5.5`, overflowing past step 99999 → `lnd*****` and a NetCDF create/permission
crash); widened to `I0.5` with a length-32 buffer.

### Scope of the restarted state (carbon, canopy, crop, lake, TOPMODEL)

The restart persists the **full internal Noah-MP prognostic state**, not just soil/snow:

- canopy / vegetation: `CANICEXY`, `CANLIQXY`, `LAI`, `XSAI`, `TVXY`, `TGXY`, `EAHXY`, `TAHXY`;
- dynamic-veg carbon pools: `LFMASSXY`, `RTMASSXY`, `STMASSXY`, `WOODXY`, `GRAINXY`, `GDDXY`,
  and soil-carbon pools `STBLCPXY`, `FASTCPXY`;
- crop (`PGSXY`), wetland (`WSURFXY`), lake (`WSLAKEXY`), TOPMODEL saturated fraction (`FSATXY`);
- the water/energy balance accumulators (`ACC_*XY`) and the integer `ITIMESTEP`.

Every one of these is `allocated()`-guarded, so it is only written/read when the scheme that
owns it is active (e.g. the carbon pools exist only under dynamic vegetation / crop). This is
what makes a restart bitwise-consistent with a continuous run for **any** namelist option set.

**Important scoping note (to avoid overclaiming):** carbon and canopy are *internal Noah-MP
prognostics that we made restart-safe*. They are **not** part of the per-step atmosphere↔land
exchange with ERF's dynamical core — ERF does not carry a CO₂/biomass tracer and Noah-MP does
not return a carbon (NEE/GPP/NPP) flux to the dycore. The per-step two-way coupling is the
surface energy/moisture/momentum/precip exchange (§1, §3); carbon/canopy state simply
evolves inside Noah-MP and is checkpointed so it survives restarts.

---

## 5. Files changed in `drivers/erf/`

- `NoahmpIOVarType.F90` — `RAINBL/SR/MP_RAINNC/MP_SNOW/MP_GRAUP` → `real(kind=C_DOUBLE)`.
- `NoahmpIO.H` — `double*` members + `NoahArray2D<double>` for the precip fields.
- `NoahmpIO.cpp` — bind each precip field `NoahArray2D<double>(fptr.X, …)`.
- `NoahmpIO_fi.F90` — `C_LOC` wiring + `type(C_PTR)` entries for the precip fields.
- `NoahmpDriverMainMod.F90` — drop the lumping remap; keep caller-supplied `MP_*`; zero only
  the absent channels; `SNOWBL=0` (inert); RAINBL still set.
- `NoahmpWriteRestartMod.F90`, `NoahmpReadRestartMod.F90` — full-state checkpoint I/O.
- `NoahmpWriteLandMod.F90` — filename width fix.

The shared core (`src/…`, `ForcingVarInTransferMod.F90`) is **unmodified** — it already
consumes `MP_*`/`SR` exactly as WRF does.

---

## 6. Validation summary

- Water/energy budget: `BalanceWaterCheck`/`BalanceEnergyCheck` pass with 0 aborts once the
  kind bug (§2) is fixed; soil recharges under precip (previously monotonic drying).
- `opt_snf=4` behavioral test (cold synthetic snowy column): with the real breakdown,
  `FrozenPrecipFrac > 0` and frozen precip reaches the core (before the fix it was 0 → all
  rain).
- `opt_snf=1` regression: identical to the prior validated run (ignores `MP_*`), confirming
  the change is safe for in-flight runs.
- Restart: continuous vs restarted trajectories match to checkpoint precision.
