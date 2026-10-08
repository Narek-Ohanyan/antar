"""Drive file ids and year ranges for MNEME's panel, shared by the extraction and the fit."""

VITALITY_TILE_SIZE_PX = 4864
VITALITY_TILE_IDS = {
    (0, 0): "1X5VF3X8k9DvJzpxVbO6Ccjfhq9tSvHLW", (0, 4864): "1iC4l1_UoG-2ZpAKP90CeM3e-1wbpZUYA",
    (4864, 0): "1CVhnSKyGkMuLQUgA4mKNoLNPvzZypZvj", (4864, 4864): "1JyNnH34ZEGSYcB4qMrthliOd35ThOSEP",
}
# Corrected export of 2026-10-07. The first file (id 1-m_xKLBqBZYY8EFhFb9C_9iCzcJcB-jv, 2026-09-24) read "disturbed" at every pixel without a recorded loss: Hansen's lossyear band is
# masked there, not zero, so the layer was masked and then exported as 0. It must not be used (see antar.io.gee_export.export_disturbance_ancillary and IMPLEMENTATION_LOG.md).
DISTURBANCE_FILE_ID = "1xYI2U_jQywz6DhdDLZn4jfvyuhs76m_7"
VITALITY_YEARS = list(range(2000, 2025))   # the full record, needed for the 10-year trailing window and the 2-year recovery window
# Years with both a determinable label (2000 + 10 <= year <= 2024 - 2) and a real climate covariate (CHELSA-daily precipitation ends 2019-12-31).
PANEL_YEARS = list(range(2010, 2020))
