# Cost and safety data

| File | Content | Source |
|---|---|---|
| `cost_ld50_donkers2017.csv` | Bulk price (EUR/kg of stable hydrate, EUR/MJ stored) and oral LD50 (rat) for 9 salt hydrates | Donkers et al. 2017, *Applied Energy* 199, 45-68, Table 3 p.52 (open access) [R13]. Prices are 2017 market prices from a supplier (Caldic), as stated in the paper, Section 6.7 |
| `ghs_pubchem.csv` | GHS signal word, pictograms and hazard statements for 48 compounds | PubChem PUG-View "GHS Classification" (aggregated ECHA C&L notifications), retrieved 11 Nov 2026. Only H-statements reported by at least 50 % of notifiers are kept |
| `material_compound_map.csv` | Which compound is used for the hazard of each database entry (storage salt; working fluid for NH3/methanol/ethanol pairs) | – |

An empty GHS row means PubChem lists no GHS classification for that compound record; it does **not** prove the material is harmless (e.g. cobalt chloride: check the supplier SDS).
Prices are only available for the salts in Donkers 2017; the price of other materials is not sourced yet.
