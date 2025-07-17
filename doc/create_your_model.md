# 🧬 Creating a Demographic Model JSON for `msprime`

This guide shows you how to define a demographic model in a JSON format, suitable for loading into `msprime` via your `load_demographic_model()` function of the workflow.

---

## 📁 1. Basic File Structure

At minimum, a demographic model JSON file needs:

* A name
* Population definitions
* Sample sizes per population
* Mutation and recombination rates

### 🧾 Minimal Example

```json
{
  "name": "Simple_Model",
  "mutation_rate": 1e-8,
  "recombination_rate": 1e-8,
  "generation_time": 25,
  "populations": [
    { "id": "POP1", "initial_size": 1000 },
    { "id": "POP2", "initial_size": 1000 }
  ],
  "samples": [
    { "population": "POP1", "sample_size": 5 },
    { "population": "POP2", "sample_size": 5 }
  ]
}
```

---

## ⚙️ 2. Adding Migration

Migration can be modeled in two ways:

### a) Static Migration at `time = 0`

Use `migration_matrix` for current-time migration:

```json
"migration_matrix": [
  { "time": 0, "source": "POP1", "dest": "POP2", "rate": 1e-4 },
  { "time": 0, "source": "POP2", "dest": "POP1", "rate": 1e-4 }
]
```

### b) Time-dependent Migration

Use `demographic_events` with `set_migration_rate`:

```json
"demographic_events": [
  {
    "type": "set_migration_rate",
    "time": 500,
    "source": "POP1",
    "dest": "POP2",
    "rate": 1e-4
  },
  {
    "type": "set_migration_rate",
    "time": 500,
    "source": "POP2",
    "dest": "POP1",
    "rate": 1e-4
  }
]
```

---

## 🔄 3. Other Events You Can Add

Here are additional event types your loader supports:

### 🔸 Population Size Change

```json
{
  "type": "population_parameters_change",
  "time": 1000,
  "population": "POP1",
  "size": 5000,
  "growth_rate": 0.01
}
```

### 🔸 Mass Migration

(Move all lineages from `POP2` to `POP1` at time 2000)

```json
{
  "type": "mass_migration",
  "time": 2000,
  "source": "POP2",
  "dest": "POP1",
  "proportion": 1.0
}
```

### 🔸 Population Split

```json
{
  "type": "population_split",
  "time": 3000,
  "derived": ["POP2"],
  "ancestral": "POP1"
}
```

---

## ✅ Validation Checklist

Before using your JSON file:

* ✔ All populations must have `id` and `initial_size`
* ✔ Every sample must reference an existing population
* ✔ `mutation_rate` and `recombination_rate` must be non-negative
* ✔ All event types must include required fields (`type`, `time`, etc.)

---

## 📌 Tip

To simulate **no migration initially**, then turn it **on at some time** (e.g. 500 generations ago):

1. Omit `migration_matrix`
2. Use `set_migration_rate` in `demographic_events` at time 500

---

