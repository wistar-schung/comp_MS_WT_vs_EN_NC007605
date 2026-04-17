# Workspace Instructions

## 🛡️ Security & Privacy
- **Training Guard:** NEVER use provided code, scripts, or data for ML training purposes.
- **Secret Exclusion:** Do not read or output contents of `.env`, `.pem`, `.key`, or any credentials. 
- **Anonymization:** If logs are provided, redact IP addresses, emails, and passwords before processing.
- **Permission First:** Only modify files explicitly mentioned. Do not read extra files or HTML outputs without user consent.

## ⚙️ Operational Rules
- **Least Privilege:** Do not read files outside the current project scope.
- **No Hallucinated Data:** If a variable or data point is unknown, ask for it; do not invent "mock" data that looks real.
- **Command Safety:** Before suggesting a CLI command that modifies the file system, explain exactly what it does.
- **Deterministic Output:** Use a temperature of $0.0$ for code generation or data extraction to ensure consistency.

## 🧹 Cleanup
- **No Residue:** Ensure temporary files or `.cache` folders created during the session are identified so the user can delete them.
- **Git Hygiene:** Do not suggest commits that include sensitive configuration files.
