#!/usr/bin/env xonsh

import click


@click.command()
@click.argument("samplefile")
@click.argument("refgenome")
@clikck.argument("fastqdir")
