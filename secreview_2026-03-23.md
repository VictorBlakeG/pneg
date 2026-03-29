# Security Review — pn-junc-sim
**Date:** 2026-03-23
**Reviewer:** Automated (Claude)
**Files Reviewed:** pn-junction-simulator.html, Dockerfile, nginx.conf

## Methodology
The following 13 categories were checked:
1. XSS vulnerabilities (innerHTML, document.write, eval, Function(), setTimeout/setInterval with strings)
2. External resource loading (CDNs/scripts without Subresource Integrity hashes)
3. Missing Content Security Policy
4. URL parameter handling without sanitization (location.search, URLSearchParams, hash params)
5. localStorage/sessionStorage usage and data exposure
6. postMessage handling without origin validation
7. Open redirects
8. Mixed content (HTTP resources on HTTPS page)
9. Sensitive data exposure (API keys, tokens, credentials)
10. Nginx configuration issues (missing security headers, directory listing, server version exposure)
11. Dockerfile issues (running as root, unnecessary packages)
12. Inline event handlers
13. Use of dangerous DOM APIs

## Findings

### External Resources Without Subresource Integrity (SRI)
**Severity:** Medium
**Location:** pn-junction-simulator.html:7-8
**Description:** Two external CDN scripts are loaded without `integrity` or `crossorigin` attributes:
- `https://cdn.jsdelivr.net/npm/chart.js@4.4.0/dist/chart.umd.min.js`
- `https://cdn.jsdelivr.net/npm/chartjs-plugin-annotation@latest/dist/chartjs-plugin-annotation.min.js`

If the CDN is compromised or a DNS hijack occurs, arbitrary JavaScript could be injected. The `@latest` tag on the annotation plugin is especially risky as it auto-resolves to whatever the current version is, meaning a malicious publish would immediately affect this app.

### Missing Content Security Policy
**Severity:** Medium
**Location:** pn-junction-simulator.html (no `<meta>` CSP tag), nginx.conf (no `Content-Security-Policy` header)
**Description:** No Content Security Policy is defined anywhere. A CSP would mitigate the impact of any XSS by restricting script sources, inline script execution, and other resource origins.

### innerHTML Used with Template Literals
**Severity:** Low
**Location:** pn-junction-simulator.html:584-585, 919-921, 926-933, 936-947, 951-954, 1442
**Description:** Multiple uses of `.innerHTML` with template literal strings. The values inserted come from internal calculation results (numeric `.toFixed()` / `.toExponential()` calls) and hardcoded material property names, not from user-controlled string input. The risk is low because the interpolated values are numeric outputs of `Math` functions, not arbitrary user text. However, using `innerHTML` at all creates a maintenance risk — future changes that interpolate string values could introduce XSS.

### URL Parameter Handling Without Sanitization
**Severity:** Low
**Location:** pn-junction-simulator.html:1676-1679
**Description:** The `theme` URL parameter is read via `URLSearchParams` and stored directly into `localStorage`. The value is then compared against the string `'dark'` to decide whether to click the theme button. Since the comparison is strict equality (`==='dark'`), only the exact string `'dark'` triggers the action. However, any arbitrary string passed via `?theme=...` is stored persistently in localStorage without validation. This is a minor data pollution issue, not exploitable for XSS.

### localStorage Usage
**Severity:** Low
**Location:** pn-junction-simulator.html:1555-1559, 1670, 1677-1678
**Description:** Two localStorage keys are used:
- `pnJunctionState`: stores simulator UI state (material, temperature, doping values, checkbox states). All values are numeric or boolean. The `loadState()` function at line 1558-1623 wraps restoration in a try/catch and sets values via `.value` and `.checked` properties (not innerHTML), which is safe.
- `isd-theme`: stores `'dark'` or `'light'`.

No sensitive data is stored. The data is all simulator configuration. Risk is minimal.

### Missing Nginx Security Headers
**Severity:** Medium
**Location:** nginx.conf
**Description:** The nginx configuration is minimal (10 lines) and lacks all recommended security headers:
- `X-Content-Type-Options: nosniff`
- `X-Frame-Options: DENY` (or `SAMEORIGIN`)
- `X-XSS-Protection: 0` (legacy, but still commonly set)
- `Referrer-Policy`
- `Permissions-Policy`
- `Strict-Transport-Security` (HSTS, if served over HTTPS)
- `Content-Security-Policy`

The `server_tokens` directive is not set to `off`, so nginx will expose its version in response headers (`Server: nginx/x.x.x`).

### Nginx Catch-All SPA Routing
**Severity:** Info
**Location:** nginx.conf:8
**Description:** `try_files $uri $uri/ /pn-junction-simulator.html` serves the HTML file for any unmatched path. This is standard SPA behavior but means that any path on this host returns the app HTML rather than a 404. This is not a vulnerability per se, but could confuse security scanners or make path-based access controls harder to reason about.

### Dockerfile Runs as Root
**Severity:** Low
**Location:** Dockerfile
**Description:** The Dockerfile does not include a `USER` directive to drop privileges. The nginx process inside the container runs as root by default. While `nginx:alpine` typically has nginx worker processes drop to a `nginx` user, the master process remains root. Best practice is to add `USER nginx` or use a non-root base image.

### Dockerfile Missing Health Check
**Severity:** Info
**Location:** Dockerfile
**Description:** No `HEALTHCHECK` instruction is defined. While not a direct security issue, health checks support container orchestration reliability.

### No postMessage, Open Redirects, Cookies, or Sensitive Data
**Severity:** Info
**Location:** N/A
**Description:** No `postMessage` listeners, no `window.open` or redirect patterns, no cookies, no API keys, tokens, or credentials were found. The app is a purely client-side physics simulator with no backend communication.

### No Mixed Content
**Severity:** Info
**Location:** N/A
**Description:** All external resources use HTTPS URLs. No HTTP resources were found.

## Recommended Fixes

A. **Add SRI hashes to CDN scripts** (pn-junction-simulator.html:7-8). Generate integrity hashes with `shasum` and add `integrity="sha384-..."` and `crossorigin="anonymous"` attributes to both `<script>` tags. Pin the annotation plugin to a specific version instead of using `@latest`.

B. **Add Content Security Policy** either as a `<meta>` tag in the HTML or (preferably) as an nginx header. A reasonable starting policy:
```
Content-Security-Policy: default-src 'self'; script-src 'self' https://cdn.jsdelivr.net 'unsafe-inline'; style-src 'self' 'unsafe-inline'; img-src 'self' data:; connect-src 'none';
```

C. **Add security headers to nginx.conf** inside the `server` block:
```nginx
add_header X-Content-Type-Options "nosniff" always;
add_header X-Frame-Options "DENY" always;
add_header Referrer-Policy "strict-origin-when-cross-origin" always;
server_tokens off;
```

D. **Add a non-root USER directive to the Dockerfile** after the COPY instructions:
```dockerfile
USER nginx
```

E. **Validate the `theme` URL parameter** before storing it in localStorage (pn-junction-simulator.html:1676-1677). Only store the value if it equals `'dark'` or `'light'`:
```javascript
var _tp = new URLSearchParams(window.location.search).get('theme');
if (_tp === 'dark' || _tp === 'light') localStorage.setItem('isd-theme', _tp);
```

F. **Consider replacing innerHTML with textContent or DOM APIs** where possible (lines 584-585, 919-954, 1442) to eliminate the class of risk entirely. Where HTML structure is needed (e.g., the calculations section), use a templating approach that escapes values.

## Summary
This is a static, client-side physics simulator with no backend API, no authentication, and no sensitive data handling. The overall risk profile is low. The most actionable findings are: (1) external CDN scripts loaded without SRI hashes, which is the primary supply-chain risk; (2) missing nginx security headers, which is a hardening gap; and (3) the Dockerfile running as root, which is a container best-practice violation. No high-severity vulnerabilities were identified.
