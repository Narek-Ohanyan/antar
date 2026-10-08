# Observed forest dieback in Armenia: what ANTAR could use, and who may hold it

MNEME (the statistical dieback hazard) is built and tested but cannot be fitted: the satellite record yields 162 dieback onsets in 570,175 forest pixel-years but cannot tell them from changes in the record itself (92 pixels with a persistent decline against 1,370 with a persistent rise) (see the Models page and `IMPLEMENTATION_LOG.md`, 2026-10-07/08). No georeferenced record of tree mortality or dieback in Armenia was found in open sources. Dated, located records of the kinds below would let the model be fitted and checked against something other than a satellite index.

## What would help, in order of value

1. **Located dieback or mortality observations with a date**: stand polygons or plot coordinates, the year (or season) the dieback was first noticed, species, and, if recorded, the cause (drought, bark beetle, nematode, fire, felling). Any year since 2000; the satellite record is 2000-2024.
2. **Repeated forest inventory plots**: tree status (alive, dead, felled) at two or more visits, with coordinates, species and diameter. The National Forest Inventory now being set up (Open Foris Arena) will provide this from its second cycle.
3. **Pest and disease surveys with coordinates**: for example the pine dieback around Lake Sevan and the pest-affected stands treated in 2021 (the Ministry of Environment gave about 12,000 ha).
4. **Tree-ring series**: ring widths of juniper and oak from Armenian sites (Opala-Owczarek et al. 2021 report chronologies of up to 140 years), or of pine and beech. They give an annual growth response to drought that can be compared with the model.
5. **Stand age or structure**: forest management inventories (Hayantar) with stand age, height and canopy cover.

## Who may hold it

| Holder | What they may have |
|---|---|
| Hayantar (State Forestry Non-Commercial Organisation) and the National Forest Inventory team (FAO / Green Climate Fund) | forest management inventories, pest records, the new inventory plots |
| Sevan National Park | the dry pine stand records since 2017 (8,269 dry trees counted; bark beetles, a suspected nematode) |
| Scientific Center of Zoology and Hydroecology, National Academy of Sciences | pest identification and survey records |
| Ministry of Environment, and the Hydrometeorology and Monitoring Center | the 2021 pest-treatment areas; monitoring of the Sevan pine dieback |
| FAO project TCP/ARM/4005 (forest pests and diseases) | surveys under the pest and disease project |
| WSL and the FORACCA team (Rueetschi, Ginzler; the remote-sensing group) | the vegetation height model, and possibly disturbance or health layers; contact through the project |
| Authors of the Armenian dendrochronology (Opala-Owczarek and colleagues, University of Silesia) | the juniper and oak chronologies |

## Short request text

> We are building a climate-resilience model for Armenian forests (ANTAR) and would like to test it against observed dieback. Could you share, under whatever terms you prefer, any located and dated records of tree dieback, mortality or pest damage in Armenian forests (since 2000), or repeated inventory plots with tree status? Coordinates (or stand polygons), the year of observation and the species are the minimum we need; the cause, if known, helps a great deal. We would credit the source in the model's documentation and share the result with you.

## What the data would be used for

- Replace the satellite-only dieback label of MNEME with a located, dated one, and estimate how often the satellite rule agrees with it (sensitivity and specificity), which also allows the misclassification correction (Rogan and Gladen 1978) that is currently not applied.
- Fit the hazard model on real events, and compare its predictions with XYLEM's mechanistic hazard.
- Check, independently, the ranking of drought stress that XYLEM produces.
