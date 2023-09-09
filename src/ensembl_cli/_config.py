import configparser
import os
import pathlib

from dataclasses import dataclass
from typing import Iterable

from ensembl_cli.species import Species, species_from_ensembl_tree


@dataclass
class Config:
    host: str
    remote_path: str
    release: str
    staging_path: os.PathLike
    install_path: os.PathLike
    species_dbs: Iterable[str]
    align_names: Iterable[str]
    tree_names: Iterable[str]

    @property
    def db_names(self) -> Iterable[str]:
        from ensembl_cli.species import Species

        for species in self.species_dbs:
            yield Species.get_ensembl_db_prefix(species)

    @property
    def staging_genomes(self):
        return self.staging_path / "genomes"

    @property
    def install_genomes(self):
        return self.install_path / "genomes"

    @property
    def staging_homologies(self):
        return self.staging_path / "compara" / "homologies"

    @property
    def install_homologies(self):
        return self.install_path / "compara" / "homologies"

    @property
    def staging_aligns(self):
        return self.staging_path / "compara" / "aligns"

    @property
    def install_aligns(self):
        return self.install_path / "compara" / "aligns"


def read_config(config_path) -> Config:
    """returns ensembl release, local path, and db specifics from the provided
    config path"""
    from ensembl_cli.download import download_ensembl_tree

    parser = configparser.ConfigParser()

    with config_path.expanduser().open() as f:
        parser.read_file(f)

    release = parser.get("release", "release")
    host = parser.get("remote path", "host")
    remote_path = parser.get("remote path", "path")
    remote_path = remote_path[:-1] if remote_path.endswith("/") else remote_path
    staging_path = (
        pathlib.Path(parser.get("local path", "staging_path")).expanduser().absolute()
    )
    install_path = (
        pathlib.Path(parser.get("local path", "install_path")).expanduser().absolute()
    )

    species_dbs = {}
    get_option = parser.get
    align_names = []
    tree_names = []
    for section in parser.sections():
        if section in ("release", "remote path", "local path"):
            continue

        if section == "compara":
            value = get_option(section, "align_names", fallback=None)
            align_names = [] if value is None else [n.strip() for n in value.split(",")]
            value = get_option(section, "tree_names", fallback=None)
            tree_names = [] if value is None else [n.strip() for n in value.split(",")]
            continue

        dbs = [db.strip() for db in get_option(section, "db").split(",")]

        # handle synonyms
        species = Species.get_species_name(section, level="raise")
        species_dbs[species] = dbs

    if tree_names:
        # add all species in the tree to species_dbs
        for tree_name in tree_names:
            tree = download_ensembl_tree(host, remote_path, release, tree_name)
            sp = species_from_ensembl_tree(tree)
            species_dbs.update(sp)

    return Config(
        host=host,
        remote_path=remote_path,
        release=release,
        staging_path=staging_path,
        install_path=install_path,
        species_dbs=species_dbs,
        align_names=align_names,
        tree_names=tree_names,
    )