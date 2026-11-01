-- Thermochemical Energy Storage (TCES) database schema (SQLite)
PRAGMA foreign_keys = ON;

CREATE TABLE IF NOT EXISTS classification (
    mechanism      TEXT PRIMARY KEY,      -- e.g. 'Salt hydrate dehydration'
    family         TEXT NOT NULL,         -- 'Sorption' | 'Chemical reaction'
    subfamily      TEXT NOT NULL          -- 'Adsorption' | 'Absorption' | 'Gas-solid' | 'Gas-gas'
);

CREATE TABLE IF NOT EXISTS references_ (
    ref_id         TEXT PRIMARY KEY,
    authors        TEXT,
    year           INTEGER,
    title          TEXT,
    journal        TEXT,
    volume_pages   TEXT,
    doi            TEXT,
    url            TEXT
);

CREATE TABLE IF NOT EXISTS species (
    species             TEXT PRIMARY KEY,
    phase               TEXT CHECK (phase IN ('s','l','g')),
    dHf_kJ_mol          REAL,
    S_J_molK            REAL,
    molar_mass_g_mol    REAL,
    density_g_cm3       REAL,
    source              TEXT,
    source_url          TEXT,
    ref_id              TEXT REFERENCES references_(ref_id),
    dHf_NBS_kJ_mol      REAL,
    S_NBS_J_molK        REAL,
    density_source      TEXT
);

CREATE TABLE IF NOT EXISTS reactions (
    id                  TEXT PRIMARY KEY,
    name                TEXT NOT NULL,
    reaction            TEXT NOT NULL,
    mechanism           TEXT REFERENCES classification(mechanism),
    working_gas         TEXT,
    data_mode           TEXT CHECK (data_mode IN ('calc','lit')),
    dH_kJ_mol_rxn       REAL,
    dS_J_molK_rxn       REAL,
    nu_gas              REAL,
    dH_kJ_per_mol_gas   REAL,
    T_eq_C              REAL,
    T_operating_C       REAL,
    T_basis             TEXT,
    temperature_class   TEXT,
    M_left_g_mol        REAL,
    density_left_g_cm3  REAL,
    ED_grav_MJ_kg       REAL,
    ED_grav_kWh_kg      REAL,
    ED_vol_GJ_m3        REAL,
    ED_vol_basis        TEXT,
    ED_bed_GJ_m3        REAL,
    T_release_C         REAL,
    volume_decrease_pct REAL,
    DRH_pct             REAL,
    cyclability         TEXT,
    main_issue          TEXT,
    verification        TEXT,
    source_location     TEXT,
    notes               TEXT
);

CREATE TABLE IF NOT EXISTS reaction_references (
    reaction_id   TEXT REFERENCES reactions(id),
    ref_id        TEXT REFERENCES references_(ref_id),
    PRIMARY KEY (reaction_id, ref_id)
);

-- Handy view: one row per reaction with its classification and DOIs
CREATE VIEW IF NOT EXISTS v_reactions_full AS
SELECT r.*, c.family, c.subfamily,
       GROUP_CONCAT(ref.doi, '; ') AS dois
FROM reactions r
JOIN classification c ON c.mechanism = r.mechanism
LEFT JOIN reaction_references rr ON rr.reaction_id = r.id
LEFT JOIN references_ ref ON ref.ref_id = rr.ref_id
GROUP BY r.id;
