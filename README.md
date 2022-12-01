# msions

'A python package for creating MS TIC and ion plots'

## Installation

```bash
$ pip install msions
```

## Usage
`msions` can be used to work with Hardklor and Kronik files and to create MS TIC and ion plots.

```python
from msions.hardklor import hk2df
from msions.hardklor import summarize_df

hk_file = "test.hk" # path to your file
hk_df = hk2df(hk_file)
sum_hk_df = summarize_df(hk_df)
```

## Contributing

Interested in contributing? Check out the contributing guidelines. Please note that this project is released with a Code of Conduct. By contributing to this project, you agree to abide by its terms.

## License

`msions` was created by Danielle Faivre. It is licensed under the terms of the Apache License 2.0 license.

## Credits

`msions` was created with [`cookiecutter`](https://cookiecutter.readthedocs.io/en/latest/) and the `py-pkgs-cookiecutter` [template](https://github.com/py-pkgs/py-pkgs-cookiecutter).
