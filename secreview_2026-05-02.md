---
name: Security Fixes Applied — pneg/ (deployed pn-junc-sim service)
description: Addendum to 2026-03-23 secreview applying the recommended fixes to the deployed copy
type: secreview
---

# Security Fixes Applied — pneg/
**Date:** 2026-05-02
**Reviewer:** Automated (Claude)
**Reference:** `secreview_2026-03-23.md` (in this directory)

This addendum applies the previously-documented recommended fixes to the deployed pn-junc-sim service files. The 2026-03-23 review remains the authoritative findings record.

## Files Modified
- `pn-junction-simulator.html`
  - Added CSP meta: `default-src 'self'; script-src 'self' https://cdn.jsdelivr.net 'unsafe-inline'; style-src 'self' 'unsafe-inline'; img-src 'self' data:; connect-src 'none'; object-src 'none'; base-uri 'self'; frame-ancestors 'none'`
  - Added `integrity="sha384-e6nUZLBkQ86NJ6TVVKAeSaK8jWa3NhkYWZFomE39AvDbQWeie9PlQqM3pmYW5d1g"` + `crossorigin="anonymous"` on chart.js@4.4.0
  - Pinned `chartjs-plugin-annotation` from `@latest` to `@3.0.1` and added `integrity="sha384-oNtu+d18330MVFpltUTve1DatxCkkctlpA2AC3GulbVFOSqhHdDat3qHse/Lbuek"` + `crossorigin="anonymous"`

- `nginx.conf`
  - Added `server_tokens off`
  - Added headers: `X-Content-Type-Options: nosniff`, `X-Frame-Options: DENY`, `Referrer-Policy: strict-origin-when-cross-origin`, `Permissions-Policy: geolocation=(), microphone=(), camera=()`

## Findings Not Addressed
- **innerHTML with internal numeric data** (lines 584-585, 919-921, 926-933, 936-947, 951-954, 1442) — Low risk, all interpolated values are numeric `.toFixed()`/`.toExponential()` outputs from physics calculations. Not refactored; covered by new CSP.
- **URL parameter `?theme` not validated before localStorage write** (lines 1676-1679) — Low risk (value is only string-compared, never rendered to DOM). Not addressed today; can be hardened later.
- **Dockerfile `USER nginx`** — Not applied. Would require chowning `/var/cache/nginx`, `/var/log/nginx`, and pid file path, OR switching to `nginxinc/nginx-unprivileged:alpine`. Should be tested in CI before applying to a deployed image.

## Other files in this directory not reviewed today
The 2026-03-23 review covered only `pn-junction-simulator.html`. This directory also contains `pn-junc-sim.html`, `pneg.html`, `pneg.jsx`, and `Mobility Calculation for PN Diodes (HTML version).html`. These were not in the original review scope and remain unreviewed.
