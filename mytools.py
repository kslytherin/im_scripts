import numpy as np
import re


def find_dict_nan(some_dict):
    nan_indices = {}
    for key_name in some_dict.keys():
        try:
            nan_TF_values = np.isnan(some_dict[key_name].values)
            if any(nan_TF_values):
                nan_i_set = np.where(nan_TF_values)
                print(
                    "{0} has {1}/{2} NaNs".format(
                        key_name, len(nan_i_set[0]), len(nan_TF_values)
                    )
                )
                nan_indices[key_name] = nan_i_set[0]
        except TypeError:
            print("TypeError issue for " + key_name)
    return nan_indices


def find_nparray_nan(some_nparray):
    nan_TF_values = np.isnan(some_nparray)
    nan_i_set = np.where(nan_TF_values)
    print("Has {0}/{1} NaNs".format(len(nan_i_set[0]), len(nan_TF_values)))
    nan_indices = nan_i_set[0]
    return nan_indices


def get_dict_key_len(some_dict):
    key_len = {}
    for key_name in some_dict.keys():
        print(key_name)
        if isinstance(some_dict[key_name], (float, bool)):
            key_len[key_name] = 1
            print("1")
            continue
        else:
            key_len[key_name] = len(some_dict[key_name])
        print(len(some_dict[key_name]))
    return key_len


def split_num_letters_caps(test_str):
    # Note: doesnt do a good job with smth like this: "AFDS1435r43"
    num_letter_tup = split_num_letters(test_str)
    split_words = split_cap_words_crop(num_letter_tup[0])
    split_words.append(num_letter_tup[1])
    return split_words


def split_num_letters(test_str):
    # Note: doesnt do a good job with smth like this: "AFDS1435r43"
    temp = re.compile("([a-zA-Z]+)([0-9]+)")
    return list(temp.match(test_str).groups())


def split_cap_words_crop(test_str):
    return re.findall("[A-Z][a-z]*", test_str)


def split_cap_words(test_str):
    return re.findall("[A-Z][^A-Z]*", test_str)
