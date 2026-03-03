"""
Distance-method hierarchy and email alert template settings.

Used by: progenitors.util
Location: progenitors/settings/util.py
"""
# Distance method preference order: name -> list of keys to match in NED/refs
DM_HIERARCHY = [
    {"name": "Cepheids", "keys": ["Cepheids", "Cepheid"]},
    {"name": "TRGB", "keys": ["TRGB"]},
    {"name": "Tully-Fisher", "keys": ["Tully-Fisher", "Tully est"]},
    {"name": "Fundamental Plane", "keys": ["FP"]},
    {"name": "SNIa", "keys": ["SNIa"]},
    {"name": "SNII", "keys": ["SNII optical", "SNII radio"]},
    {"name": "Sosies", "keys": ["Sosies"]},
    {"name": "Planetary Nebula Luminosity Function", "keys": ["PNLF"]},
    {"name": "Ring Diameter", "keys": ["Ring Diameter"]},
]

# Email alert (from-addr label, SMTP, subject); to/login/password come from env
EMAIL_FROM_ADDR = "Supernova Progenitor Alerts"
EMAIL_SMTP_SERVER = "smtp.gmail.com"
EMAIL_SMTP_PORT = 587
EMAIL_SUBJECT = "Supernova Progenitor Target Summary"
EMAIL_MSG_TEMPLATE = """<html><body><p>Bright transients for APF follow up</p>
{targets}
<p>CDK</p>
</body></html>"""
