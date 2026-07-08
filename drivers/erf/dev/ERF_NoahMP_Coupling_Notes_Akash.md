# ERF–Noah-MP Coupling: Design Notes for the ERF Driver

**Scope.** This document describes the ERF driver of Noah-MP (`drivers/erf/`) and the
accompanying changes on the ERF side (`Source/LandSurfaceModel/Noah-MP/ERF_NOAHMP.{H,cpp}`).
It records the coupling contract that spans the C++/Fortran boundary — precipitation forcing,
surface-flux units, the field-kind requirement, and restart completeness — together with the
physical and software rationale for each design decision. The ERF driver and the ERF C++
source must be read together, because the precipitation and flux conventions are shared
across the boundary.

**Relation to prior work.** This effort parallels erf-model/ERF PR #3405, *"Expose IO and
precipitation variables from Noah-MP, enhance ERF–Noah-MP coupling for precipitation."* The
two share the same objective and agree on the central structural decision (surface fluxes are
kinematic, i.e. divided by density). This document additionally documents three corrections
not present in #3405; these are summarized in §1 and detailed in the following sections.

---

## Contents

1. [Summary of differences from PR #3405](#1-summary-of-differences-from-pr-3405)
2. [Precipitation forcing into the land model](#2-precipitation-forcing-into-the-land-model)
3. [Field-kind consistency across the C/Fortran boundary](#3-field-kind-consistency-across-the-cfortran-boundary)
4. [Surface-flux units returned to ERF](#4-surface-flux-units-returned-to-erf)
5. [Restart and checkpoint completeness](#5-restart-and-checkpoint-completeness)
6. [Files modified](#6-files-modified)
7. [Validation](#7-validation)
8. [Appendix: evidence that SNOWBL is inactive in the WRF driver](#8-appendix-evidence-that-snowbl-is-inactive-in-the-wrf-driver)

---

## 1. Summary of differences from PR #3405

Both PR #3405 and this work feed precipitation into Noah-MP and expose its I/O across the ERF
boundary, and both adopt the kinematic (density-divided) surface-flux convention. This branch
introduces three additional corrections:

| Aspect | PR #3405 | This work | Consequence |
|---|---|---|---|
| Momentum and moisture flux, division by ρ | `HFX/(ρ·Cp)`, `LH/(ρ·L)`, `τ/ρ` | Identical | Consistent |
| Exner factor on the heat flux | Absent; supplies `HFX/(ρ·Cp)`, a kinematic *temperature* flux | Present; supplies `HFX/(ρ·Cp)·(p₀/p)^(Rd/Cp)`, a kinematic *potential-temperature* flux | ERF advances potential temperature; omitting Exner biases surface heating |
| `opt_snf = 4` rain/snow partition | Lumped: `MP_SNOW = MP_GRAUP = 0` in the driver | Microphysics `MP_RAINNC`, `MP_SNOW`, `MP_GRAUP` supplied by the caller | With the lumped path, the frozen fraction is identically zero, so all precipitation is treated as rain |
| Kind of `RAINBL`, `SR`, `MP_*` | Bound to C `double*` | Declared `real(kind=C_DOUBLE)` | A 4-byte field bound to `double*` is reinterpreted and read as invalid data |
| Restart completeness | Out of scope | Full prognostic and scheme-gated state checkpoint I/O | Restarts continue rather than cold-start the land model |

The remaining sections provide the detail behind this table.

---

## 2. Precipitation forcing into the land model

### 2.1 Prior behavior

The ERF driver did not supply precipitation to Noah-MP: `RAINBL` was never set, so the land
model integrated with zero precipitation and soil moisture could only decrease over an event.
The microphysics breakdown fields were derived inside `NoahmpDriverMainMod.F90` by lumping:

```fortran
RAINNCV    = RAINBL          ! RAINBL was effectively zero
SNOWNCV    = SNOWBL          ! SNOWBL uninitialized/zero
GRAUPELNCV = 0.0
...
MP_RAINNC  = RAINNCV
MP_SNOW    = SNOWNCV         ! -> 0
MP_GRAUP   = GRAUPELNCV      ! -> 0
```

PR #3405 retains this lumping (`MP_SNOW = SNOWNCV`, `MP_GRAUP = GRAUPELNCV`, both zero).

### 2.2 Reference behavior in WRF

WRF's Noah-MP caller (`drivers/wrf/NoahmpWRFmainMod.F90`, lines 560–568) assigns
`MP_RAINNC`, `MP_SNOW`, `MP_GRAUP`, `MP_RAINC`, `MP_SHCV`, and `MP_HAIL` directly from the
microphysics, with snow and graupel as subsets of the total non-convective amount
`MP_RAINNC`. The Noah-MP driver does not re-derive these fields.

### 2.3 Implementation

The ERF path is aligned with WRF's caller:

- The ERF C++ driver computes, per land-model call interval, the accumulated-precipitation
  increments in millimetres from the microphysics accumulators and passes them across the
  boundary:
  - `RAINBL   = Δrain_accum`  — total surface precipitation (mm), clamped to a physical ceiling
  - `SR       = Δfrozen / Δtotal`  — frozen fraction, used by the auxiliary precipitation path
  - `MP_RAINNC = Δrain_accum`  — total non-convective precipitation (mm)
  - `MP_SNOW   = Δsnow_accum`  — snow-plus-ice subset (mm)
  - `MP_GRAUP  = Δgraup_accum` — graupel subset (mm)

  ERF's Morrison microphysics provides these with the same subset convention: `rain_accum`
  (total), `snow_accum` (ice + snow), and `graup_accum` (graupel).

- In `NoahmpDriverMainMod.F90`, the lumping remap is removed so that the caller-supplied
  `MP_*` fields are not overwritten. Only the channels ERF does not represent — `MP_RAINC`,
  `MP_SHCV`, and `MP_HAIL` (convective, shallow-convective, and hail) — are set to zero, which
  matches WRF when those optional arguments are absent. `RAINBL` (the total) is still set and
  continues to feed the `opt_snf = 1` path and `ACSNOW`.

### 2.4 Rationale: the `opt_snf = 4` partition

The shared Noah-MP core (`src/AtmosForcingMod.F90`, `OptRainSnowPartition == 4`) computes the
frozen fraction as

```
FrozenPrecipFrac = PrecipFrozenTot / PrecipNonConvRefHeight
PrecipFrozenTot  = (MP_SNOW + MP_GRAUP + MP_HAIL) / DTBL
```

If `MP_SNOW = MP_GRAUP = 0`, then `FrozenPrecipFrac = 0` and all precipitation is treated as
rain regardless of temperature. Supplying the microphysics breakdown corrects this. The
`opt_snf = 1` (temperature-based) partition does not consult `MP_*`, so this change is
behavior-preserving for existing `opt_snf = 1` configurations.

`SNOWBL` is not the mechanism by which snow reaches the core; it is an inactive field in the
WRF driver. It is set to zero here defensively, because an uninitialized `SNOWBL` propagated
into `MP_SNOW` and produced a water-balance abort (see §3). Supporting evidence is given in
the [appendix](#8-appendix-evidence-that-snowbl-is-inactive-in-the-wrf-driver).

---

## 3. Field-kind consistency across the C/Fortran boundary

`RAINBL`, `SR`, and `MP_RAINNC`/`MP_SNOW`/`MP_GRAUP` cross the C++/Fortran boundary as C
`double*` (declared in `NoahmpIO.H`, wired with `C_LOC` in `NoahmpIO_fi.F90`). In this build
`kind_noahmp` is single precision (`DOUBLE_PREC` undefined), so a field declared
`real(kind=kind_noahmp)` occupies four bytes. Binding a 4-byte array to a `double*`
reinterprets the underlying bytes, and Noah-MP reads invalid values (observed magnitudes of
10¹⁰–10²⁵), which then fail the water-balance check.

**Resolution.** Every field exposed to C++ is declared `real(kind=C_DOUBLE)` in
`NoahmpIOVarType.F90` (`RAINBL`, `SR`, `MP_RAINNC`, `MP_SNOW`, `MP_GRAUP`). This is a
prerequisite for any trustworthy precipitation value. Any field subsequently added to the C++
boundary must likewise be declared `C_DOUBLE`.

---

## 4. Surface-flux units returned to ERF

This section concerns the ERF C++ side (`ERF_NOAHMP.cpp`); it is documented here because the
convention is part of the shared coupling contract.

Noah-MP returns sensible and latent heat fluxes `HFX`, `LH` in W m⁻² and momentum fluxes
`TAU_EW`, `TAU_NS` in N m⁻² (density included). ERF advances `RhoTheta` (ρ·θ) and consumes
surface fluxes in the kinematic Monin–Obukhov convention — potential temperature, without a
density factor. This is established by direct inspection of the ERF code that consumes the
fluxes:

- `ERF_MOSTStress.H` documents `surf_temp_flux` as `⟨θ′w′⟩` (K m s⁻¹) and forms
  `t_star = −surf_temp_flux / u_star`, with no density factor.
- The PBL schemes (`ERF_ComputeDiffusivityMYNN25/MYNNEDMF/MRF`, `ERF_PBLHeightYSU_MRF.H`)
  form the Obukhov length `L = −θ·u*³ / (κ·g·⟨θ′w′⟩)`, which is dimensionally consistent only
  if the heat flux is the kinematic potential-temperature flux.
- The land-to-surface-layer bridge (`ERF_SurfaceLayer.cpp`) computes
  `u_star = sqrt(tau13² + tau23²)`, so the stress components must be `u*²` (m² s⁻²), i.e.
  `TAU/ρ`.

The conversion applied on this branch is therefore

```
t_flux = HFX/(ρ·Cp_d) · (p₀/p)^(Rd/Cp_d)     kinematic potential-temperature flux [K m s⁻¹]
q_flux = LH /(ρ·L_v)                          kinematic moisture flux [kg kg⁻¹ m s⁻¹]
tau13  = TAU_EW/ρ,   tau23 = TAU_NS/ρ         kinematic momentum flux [m² s⁻²]
```

Here `p` is the full pressure obtained from the equation of state applied to `RhoTheta`
(`getPgivenRTh`), not a perturbation; `p₀` is the reference pressure (10⁵ Pa). The Exner
factor is close to unity near the surface and increases with height.

**Difference from PR #3405.** PR #3405 applies the density division correctly but omits the
Exner factor, supplying `HFX/(ρ·Cp)` — a kinematic *temperature* flux — where ERF's energy
variable is *potential* temperature. The missing factor `(p₀/p)^(Rd/Cp)` biases the surface
heat flux (on the order of one percent near the surface, larger aloft). This branch adds the
Exner factor to the heat flux only.

This convention is required by ERF's own surface-layer and PBL discretization; it is not a
reproduction of WRF's internal convention.

*Development note.* An intermediate revision of this work incorrectly removed the density
division (reasoning from the interior diffusion operator `−ρα ∂θ/∂z` alone). Review of the
PBL and Monin–Obukhov consumers established that the kinematic form is required, so the
density division was retained; the only genuine defect in the original flux code was the
missing Exner factor. The net change relative to the original code is the addition of the
Exner factor on the heat flux.

---

## 5. Restart and checkpoint completeness

Earlier ERF checkpoints persisted only a subset of Noah-MP state, so a restart cold-started
the remainder and diverged from a continuous integration. This branch adds full
prognostic-state checkpoint I/O (`NoahmpWriteRestartMod.F90`, `NoahmpReadRestartMod.F90`),
together with the balance accumulators (`ACC_*XY`) and the integer step counter `ITIMESTEP`.
All fields are guarded by `allocated()`, so the restart is correct for any namelist option
set: scheme-gated fields are written and read only when their scheme is active. The result is
a restart that reproduces a continuous run to checkpoint precision. (Related to ERF #3255.)

A separate correction in `NoahmpWriteLandMod.F90` addresses a filename-width overflow: the
timestamp format was `I5.5`, which overflows beyond step 99999 and produces an invalid
filename and a NetCDF creation failure. It is widened to `I0.5` with a length-32 buffer.

### 5.1 Scope of the persisted state

The restart persists the full internal Noah-MP prognostic state:

- Canopy and vegetation: `CANICEXY`, `CANLIQXY`, `LAI`, `XSAI`, `TVXY`, `TGXY`, `EAHXY`,
  `TAHXY`.
- Dynamic-vegetation carbon pools: `LFMASSXY`, `RTMASSXY`, `STMASSXY`, `WOODXY`, `GRAINXY`,
  `GDDXY`; soil-carbon pools `STBLCPXY`, `FASTCPXY`.
- Crop (`PGSXY`), wetland (`WSURFXY`), lake (`WSLAKEXY`), and TOPMODEL saturated fraction
  (`FSATXY`).
- Water and energy balance accumulators (`ACC_*XY`) and the step counter `ITIMESTEP`.

Each field is `allocated()`-guarded and is therefore persisted only when the owning scheme is
active (for example, the carbon pools exist only under dynamic vegetation or crop modeling).

**Scope clarification.** The carbon, canopy, crop, lake, and TOPMODEL fields are internal
Noah-MP prognostics that this work makes restart-safe. They are not part of the per-time-step
exchange with ERF's dynamical core: ERF carries no carbon or biomass tracer, and Noah-MP does
not return a carbon flux (NEE/GPP/NPP) to the dynamical core. The per-step, two-way coupling
consists of the surface energy, moisture, momentum, and precipitation exchange described in
§2 and §4; the carbon and canopy state evolves within Noah-MP and is checkpointed so that it
survives restarts.

---

## 6. Files modified

In `drivers/erf/`:

- `NoahmpIOVarType.F90` — `RAINBL`, `SR`, `MP_RAINNC`, `MP_SNOW`, `MP_GRAUP` declared
  `real(kind=C_DOUBLE)`.
- `NoahmpIO.H` — `double*` members and `NoahArray2D<double>` accessors for the precipitation
  fields.
- `NoahmpIO.cpp` — bindings `NoahArray2D<double>(fptr.X, …)` for each precipitation field.
- `NoahmpIO_fi.F90` — `C_LOC` wiring and `type(C_PTR)` entries for the precipitation fields.
- `NoahmpDriverMainMod.F90` — removal of the lumping remap; caller-supplied `MP_*` retained;
  only the absent channels zeroed; `SNOWBL` set to zero (inactive); `RAINBL` still set.
- `NoahmpWriteRestartMod.F90`, `NoahmpReadRestartMod.F90` — full-state checkpoint I/O.
- `NoahmpWriteLandMod.F90` — filename-width correction.

The shared core (`src/…`, `ForcingVarInTransferMod.F90`) is unmodified; it already consumes
`MP_*` and `SR` in the same manner as WRF.

---

## 7. Validation

- **Water and energy budget.** `BalanceWaterCheck` and `BalanceEnergyCheck` pass with zero
  aborts once the field-kind defect (§3) is corrected; soil moisture recharges under
  precipitation, in contrast to the prior monotonic drying.
- **`opt_snf = 4` behavior.** In a cold synthetic snowy-column test, the microphysics
  breakdown yields `FrozenPrecipFrac > 0` and frozen precipitation reaches the core; prior to
  the fix the fraction was zero and all precipitation was treated as rain.
- **`opt_snf = 1` regression.** Results are unchanged relative to the previously validated
  run (this path does not consult `MP_*`), confirming the change is safe for existing
  configurations.
- **Restart.** Continuous and restarted trajectories agree to checkpoint precision.

---

## 8. Appendix: evidence that SNOWBL is inactive in the WRF driver

The declaration comment on `SNOWBL` ("snow entering land model") states an intended role, but
a declaration does not constitute a use. The field is inactive in the WRF driver: an
exhaustive search (`grep -rniE SNOWBL`) over the Noah-MP tree shows that within `drivers/wrf/`
it appears in exactly three locations, all non-executable:

| Location | Role |
|---|---|
| `drivers/wrf/NoahmpIOVarType.F90:101` | type declaration |
| `drivers/wrf/NoahmpIOVarInitMod.F90:56` | allocation |
| `drivers/wrf/NoahmpIOVarInitMod.F90:595` | initialization to `undefined_real` |

There is no assignment of a computed value and no read anywhere under `drivers/wrf/`. Two
further observations from the same search corroborate this:

1. The shared core (`src/`) contains no occurrence of `SNOWBL`; the core cannot reference it.
2. The forcing-transfer module `drivers/wrf/ForcingVarInTransferMod.F90` does not reference
   `SNOWBL`. Snow enters the core there through `MP_SNOW`, `MP_GRAUP`, and `MP_RAINNC`
   (each divided by `DTBL`), for example `PrecipSnowRefHeight = MP_SNOW(I,J)/DTBL`.

In WRF, therefore, snow reaches the core through `MP_SNOW`, and `SNOWBL` is declared and
initialized but never used. By contrast, the pre-existing ERF driver did read it
(`drivers/erf/NoahmpDriverMainMod.F90:58,63`: `SNOWBL = SNOWBL·DTBL`, `SNOWNCV = SNOWBL`),
which is the lumping remap removed by this work and the source of the water-balance abort when
`SNOWBL` was uninitialized.

This result is identical at both the canonical Noah-MP development commit and the commit that
WRF ships.

**Source references.**

- Canonical Noah-MP (`NCAR/noahmp`):
  [`drivers/wrf/NoahmpIOVarType.F90#L101`](https://github.com/NCAR/noahmp/blob/badab7b4b51710037fc87f3dbf329b6be59b1b5a/drivers/wrf/NoahmpIOVarType.F90#L101)
  (declaration),
  [`NoahmpIOVarInitMod.F90#L56`](https://github.com/NCAR/noahmp/blob/badab7b4b51710037fc87f3dbf329b6be59b1b5a/drivers/wrf/NoahmpIOVarInitMod.F90#L56)
  (allocation),
  [`#L595`](https://github.com/NCAR/noahmp/blob/badab7b4b51710037fc87f3dbf329b6be59b1b5a/drivers/wrf/NoahmpIOVarInitMod.F90#L595)
  (initialization), and the active snow path
  [`ForcingVarInTransferMod.F90#L53-L56`](https://github.com/NCAR/noahmp/blob/badab7b4b51710037fc87f3dbf329b6be59b1b5a/drivers/wrf/ForcingVarInTransferMod.F90#L53-L56).
- WRF (`wrf-model/WRF`) vendors Noah-MP as the submodule `phys/noahmp → NCAR/noahmp`
  ([`.gitmodules`](https://github.com/wrf-model/WRF/blob/master/.gitmodules)), pinned at
  commit `5da0b241`; the shipped file is
  [`phys/noahmp/drivers/wrf/NoahmpIOVarType.F90#L101`](https://github.com/NCAR/noahmp/blob/5da0b241e48ecfd9a2a1bd667ed554765856d589/drivers/wrf/NoahmpIOVarType.F90#L101)
  (declaration; allocation and initialization at
  [`NoahmpIOVarInitMod.F90#L56`](https://github.com/NCAR/noahmp/blob/5da0b241e48ecfd9a2a1bd667ed554765856d589/drivers/wrf/NoahmpIOVarInitMod.F90#L56)
  and
  [`#L595`](https://github.com/NCAR/noahmp/blob/5da0b241e48ecfd9a2a1bd667ed554765856d589/drivers/wrf/NoahmpIOVarInitMod.F90#L595);
  snow path at
  [`ForcingVarInTransferMod.F90#L53-L56`](https://github.com/NCAR/noahmp/blob/5da0b241e48ecfd9a2a1bd667ed554765856d589/drivers/wrf/ForcingVarInTransferMod.F90#L53-L56)).
