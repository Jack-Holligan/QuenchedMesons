from glob import glob
import logging

from joblib import Memory
import re
import yaml

memory = Memory("cache")


class ScrapeError(RuntimeError):
    pass


def boolify(text):
    if text == "True":
        return True
    elif text == "False":
        return False
    else:
        raise ValueError(f"{text} is neither true nor false.")


def get(line, converter):
    return converter(line.split()[2].strip('"}],'))


@memory.cache
def read_file_by_name(filename):
    with open(filename, "r") as f:
        try:
            return read_file(f)
        except ScrapeError:
            logging.warning(f"Unable to read {filename}")
            return None


def read_file(f):
    data = {}
    get_mass = False
    get_seed = False
    get_ip = False
    get_fp = False

    for line in f:
        for key, converter in [
            ("NT", int),
            ("NS", int),
            ("beta", float),
            ("Digitsb", int),
            ("DPb", int),
            ("Digitsm", int),
            ("DPm", int),
            ("Nc", int),
            ("Rep", str),
        ]:
            if line.startswith(f'  RowBox[{{"{key}", ":=",'):
                try:
                    data[key] = get(line, converter)
                except IndexError as ex:
                    raise ScrapeError(ex)
                break

        if get_mass:
            get_mass = False
            match len(split_line := line.split()):
                case 2:
                    data["m"] = -float(split_line[1].strip('"}],'))
                case 1:
                    data["m"] = float(line.split('"')[1])
                case _:
                    breakpoint()
                    raise ValueError("I can't cope with this.")
        elif line.startswith('  RowBox[{"m", ":=",'):
            match len(line.split()):
                case 3:
                    data["m"] = get(line, float)
                case 2:
                    get_mass = True
                case _:
                    breakpoint()
                    raise ValueError("I can't cope with this.")

        if get_seed:
            get_seed = False
            data["SeedRandom"] = int(line.split('"')[1])
        elif line.startswith('  RowBox[{"SeedRandom", "[",'):
            get_seed = True

        if line.strip().startswith('RowBox[{"GradientFlow", "=",'):
            data["GradientFlow"] = get(line, float)

        if (
            line.startswith('  RowBox[{"Use') and line.split()[1] == '"=",'
        ) or line.startswith('  RowBox[{"Automatic'):
            split_line = line.split('"')
            data[split_line[1]] = boolify(split_line[5])

        if (
            (
                line.startswith('       RowBox[{"IP')
                or line.startswith('       RowBox[{"FP')
            )
            and (split_line := line.split())[1] == '"=",'
            and len(split_line) > 2
        ):
            split_line = line.split('"')
            data[split_line[1]] = int(split_line[5])

        # if get_ip == "skip_one":
        #     get_ip = True
        # elif get_ip is True:
        #     breakpoint()
        #     data[ip_to_get] = int(line.split()[1].strip('",'))
        #     get_ip = False
        #     del ip_to_get
        # elif line.startswith('       RowBox[{"IP') and line.split()[1] == '"=",':
        #     get_ip = "skip_one"
        #     ip_to_get = line.split('"')[1]

        # if get_fp == "skip_one":
        #     get_fp = True
        # elif get_fp is True:
        #     data[fp_to_get] = int(line.split()[1].strip('",'))
        #     get_fp = False
        #     del fp_to_get
        # elif line.startswith('       RowBox[{"FP') and line.split()[1] == '"=",':
        #     get_fp = "skip_one"
        #     fp_to_get = line.split('"')[1]
        #
    return data


def clean_filename(filename):
    Nc, NS, NT, beta, m, rep = re.match(
        r"originals/Sp([468])_S([0-9]+)T([0-9]+)B([0-9.]+)_([-0-9/.]+)_([A-Z]+)\.nb",
        filename,
    ).groups()
    return f"Nc{Nc}_S{NS}T{NT}B{beta}_m{rep}{m}"


all_data = {
    clean_filename(filename): read_file_by_name(filename)
    for filename in glob("originals/*.nb")
}

with open("ensembles.yaml", "w") as f:
    f.write(yaml.dump(all_data, Dumper=yaml.Dumper))
