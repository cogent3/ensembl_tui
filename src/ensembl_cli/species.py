import os
import re

from collections import defaultdict

from cogent3 import load_table, open_
from cogent3.util.table import Table

from .util import ENSEMBLDBRC, CaseInsensitiveString, get_resource_path


_invalid_chars = re.compile("[^a-zA-Z _]")


def load_species(species_path):
    """returns [[latin_name, common_name],..] from species_path

    if species_path does not exist, defaults to default one"""
    if not os.path.exists(species_path):
        species_path = get_resource_path("species.tsv")

    table = load_table(species_path)
    return table.tolist()


_species_common_map = load_species(os.path.join(ENSEMBLDBRC, "species.tsv"))


class SpeciesNameMap:
    """mapping between common names and latin names"""

    def __init__(self, species_common=_species_common_map):
        """provides latin name:common name mappings"""
        self._species_common = {}
        self._common_species = {}
        self._species_ensembl = {}
        self._ensembl_species = {}
        for names in species_common:
            names = list(map(CaseInsensitiveString, names))
            self.amend_species(*names)

    def __str__(self) -> str:
        return str(self.to_table())

    def __repr__(self) -> str:
        return repr(self.to_table())

    def __contains__(self, item):
        return any(
            item in attr
            for attr in (
                self._species_common,
                self._common_species,
                self._ensembl_species,
            )
        )

    def _repr_html_(self) -> str:
        table = self.to_table()
        return table._repr_html_()

    def get_common_name(self, name: str, level="raise") -> str:
        """returns the common name for the given name (which can be either a
        species name or the ensembl version)"""
        name = CaseInsensitiveString(name)
        if name in self._ensembl_species:
            name = self._ensembl_species[name]

        if name in self._species_common:
            common_name = self._species_common[name]
        elif name in self._common_species:
            common_name = name
        else:
            common_name = None

        if common_name is None:
            msg = f"Unknown species name: {name}"
            if level == "raise":
                raise ValueError(msg)
            elif level == "warn":
                print(f"WARN: {msg}")

        return str(common_name)

    def get_species_name(self, name: str, level="ignore") -> str:
        """returns the species name for the given common name"""
        name = CaseInsensitiveString(name)
        if name in self._species_common:
            return str(name)

        species_name = None
        level = level.lower().strip()
        for data in [self._common_species, self._ensembl_species]:
            if name in data:
                species_name = data[name]
        if species_name is None:
            msg = f"Unknown common name: {name}"
            if level == "raise":
                raise ValueError(msg)
            elif level == "warn":
                print(f"WARN: {msg}")
        return str(species_name)

    def get_species_names(self):
        """returns the list of species names"""
        names = sorted(self._species_common.keys())
        return [str(n) for n in names]

    def get_ensembl_db_prefix(self, name):
        """returns a string of the species name in the format used by
        ensembl"""
        name = CaseInsensitiveString(name)
        if name in self._common_species:
            name = self._common_species[name]
        try:
            species_name = self.get_species_name(name, level="raise")
        except ValueError as e:
            if name not in self._species_common:
                raise ValueError(f"Unknown name {name}") from e
            species_name = name

        return str(species_name.lower().replace(" ", "_"))

    def _purge_species(self, species_name):
        """removes a species record"""
        species_name = CaseInsensitiveString(species_name)
        if species_name not in self._species_common:
            return
        common_name = self._species_common.pop(species_name)
        ensembl_name = self._species_ensembl.pop(species_name)
        self._ensembl_species.pop(ensembl_name)
        self._common_species.pop(common_name)

    def amend_species(self, species_name, common_name):
        """add a new species, and common name"""
        species_name = CaseInsensitiveString(species_name)
        common_name = CaseInsensitiveString(common_name)
        assert "_" not in species_name, "'_' in species_name, not a Latin name?"
        self._purge_species(species_name)  # remove if existing
        self._species_common[species_name] = common_name
        self._common_species[common_name] = species_name
        ensembl_name = species_name.lower().replace(" ", "_")
        self._species_ensembl[species_name] = ensembl_name
        self._ensembl_species[ensembl_name] = species_name

    def to_table(self):
        """returns cogent3 Table"""
        rows = []
        for common in self._common_species:
            species = self._common_species[common]
            ensembl = self._species_ensembl[species]

            rows += [[species, common, ensembl]]
        return Table(
            [
                "Species name",
                "Common name",
                "Ensembl Db Prefix",
            ],
            data=rows,
            space=2,
        ).sorted()


Species = SpeciesNameMap()
