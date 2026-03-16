#!/usr/bin/env python3
import random
import string


def rng_string_gen(n_char=8):
    use_chars = string.ascii_letters + string.digits
    return ''.join(random.choice(use_chars) for _ in range(n_char))
