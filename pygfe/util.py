from typing import GenericMeta
from datetime import datetime, date
from six import integer_types, iteritems


def _deserialize(data, klass):
    """
    Deserializes dict, list, str into an object.

    :param data: dict, list or str.
    :param klass: class literal, or string of class name.

    :return: object.
    """
    if data is None:
        return None

    if klass in integer_types or klass in (float, str, bool):
        return _deserialize_primitive(data, klass)
    elif klass == object:
        return _deserialize_object(data)
    elif klass == date:
        return deserialize_date(data)
    elif klass == datetime:
        return deserialize_datetime(data)
    elif type(klass) == GenericMeta:
        if klass.__extra__ == list:
            return _deserialize_list(data, klass.__args__[0])
        if klass.__extra__ == dict:
            return _deserialize_dict(data, klass.__args__[1])
    else:
        return deserialize_model(data, klass)


def _deserialize_primitive(data, klass):
    """
    Deserializes to primitive type.

    :param data: data to deserialize.
    :param klass: class literal.

    :return: int, long, float, str, bool.
    :rtype: int | long | float | str | bool
    """
    try:
        value = klass(data)
    except UnicodeEncodeError:
        value = unicode(data)
    except TypeError:
        value = data
    return value


def _deserialize_object(value):
    """
    Return a original value.

    :return: object.
    """
    return value


def deserialize_date(string):
    """
    Deserializes string to date.

    :param string: str.
    :type string: str
    :return: date.
    :rtype: date
    """
    try:
        from dateutil.parser import parse
        return parse(string).date()
    except ImportError:
        return string


def deserialize_datetime(string):
    """
    Deserializes string to datetime.

    The string should be in iso8601 datetime format.

    :param string: str.
    :type string: str
    :return: datetime.
    :rtype: datetime
    """
    try:
        from dateutil.parser import parse
        return parse(string)
    except ImportError:
        return string


def deserialize_model(data, klass):
    """
    Deserializes list or dict to model.

    :param data: dict, list.
    :type data: dict | list
    :param klass: class literal.
    :return: model object.
    """
    instance = klass()

    if not instance.swagger_types:
        return data

    for attr, attr_type in iteritems(instance.swagger_types):
        if data is not None \
                and instance.attribute_map[attr] in data \
                and isinstance(data, (list, dict)):
            value = data[instance.attribute_map[attr]]
            setattr(instance, attr, _deserialize(value, attr_type))

    return instance


def _deserialize_list(data, boxed_type):
    """
    Deserializes a list and its elements.

    :param data: list to deserialize.
    :type data: list
    :param boxed_type: class literal.

    :return: deserialized list.
    :rtype: list
    """
    return [_deserialize(sub_data, boxed_type)
            for sub_data in data]


def _deserialize_dict(data, boxed_type):
    """
    Deserializes a dict and its elements.

    :param data: dict to deserialize.
    :type data: dict
    :param boxed_type: class literal.

    :return: deserialized dict.
    :rtype: dict
    """
    return {k: _deserialize(v, boxed_type)
            for k, v in iteritems(data)}


def get_structures():
    return {
        'KIR2DP1':
        {
            'five_prime_UTR': 1,
            'exon-1': 2,
            'intron-1': 3,
            'exon-2': 4,
            'intron-2': 5,
            'exon-3': 6,
            'intron-3': 7,
            'exon-4': 8,
            'intron-4': 9,
            'exon-5': 10,
            'intron-5': 11,
            'exon-6': 12,
            'intron-6': 13,
            'exon-7': 14,
            'intron-7': 15,
            'exon-8': 16,
            'intron-8': 17,
            'exon-9': 18,
            'intron-9': 19,
            'three_prime_UTR': 20
        },
        'KIR2DL5A':
        {
            'five_prime_UTR': 1,
            'exon-1': 2,
            'intron-1': 3,
            'exon-2': 4,
            'intron-2': 5,
            'exon-3': 6,
            'intron-3': 7,
            'exon-4': 8,
            'intron-4': 9,
            'exon-5': 10,
            'intron-5': 11,
            'exon-6': 12,
            'intron-6': 13,
            'exon-7': 14,
            'intron-7': 15,
            'exon-8': 16,
            'intron-8': 17,
            'exon-9': 18,
            'intron-9': 19,
            'three_prime_UTR': 20
        },
        'KIR2DS4':
        {
            'five_prime_UTR': 1,
            'exon-1': 2,
            'intron-1': 3,
            'exon-2': 4,
            'intron-2': 5,
            'exon-3': 6,
            'intron-3': 7,
            'exon-4': 8,
            'intron-4': 9,
            'exon-5': 10,
            'intron-5': 11,
            'exon-6': 12,
            'intron-6': 13,
            'exon-7': 14,
            'intron-7': 15,
            'exon-8': 16,
            'intron-8': 17,
            'exon-9': 18,
            'intron-9': 19,
            'three_prime_UTR': 20
        },
        'HLA-DPB1':
        {
            'five_prime_UTR': 1,
            'exon-1': 2,
            'intron-1': 3,
            'exon-2': 4,
            'intron-2': 5,
            'exon-3': 6,
            'intron-3': 7,
            'exon-4': 8,
            'intron-4': 9,
            'exon-5': 10,
            'three_prime_UTR': 11
        },
        'KIR2DS2':
        {
            'five_prime_UTR': 1,
            'exon-1': 2,
            'intron-1': 3,
            'exon-2': 4,
            'intron-2': 5,
            'exon-3': 6,
            'intron-3': 7,
            'exon-4': 8,
            'intron-4': 9,
            'exon-5': 10,
            'intron-5': 11,
            'exon-6': 12,
            'intron-6': 13,
            'exon-7': 14,
            'intron-7': 15,
            'exon-8': 16,
            'intron-8': 17,
            'exon-9': 18,
            'intron-9': 19,
            'three_prime_UTR': 20
        },
        'KIR3DP1':
        {
            'five_prime_UTR': 1,
            'exon-1': 2,
            'intron-1': 3,
            'exon-2': 4,
            'intron-2': 5,
            'exon-3': 6,
            'intron-3': 7,
            'exon-4': 8,
            'intron-4': 9,
            'exon-5': 10,
            'intron-5': 11,
            'exon-6': 12,
            'intron-6': 13,
            'exon-7': 14,
            'intron-7': 15,
            'exon-8': 16,
            'intron-8': 17,
            'exon-9': 18,
            'intron-9': 19,
            'three_prime_UTR': 20
        },
        'HLA-DRB4':
        {
            'five_prime_UTR': 1,
            'exon-1': 2,
            'intron-1': 3,
            'exon-2': 4,
            'intron-2': 5,
            'exon-3': 6,
            'intron-3': 7,
            'exon-4': 8,
            'intron-4': 9,
            'exon-5': 10,
            'intron-5': 11,
            'exon-6': 12,
            'three_prime_UTR': 13
        },
        'KIR2DS5':
        {
            'five_prime_UTR': 1,
            'exon-1': 2,
            'intron-1': 3,
            'exon-2': 4,
            'intron-2': 5,
            'exon-3': 6,
            'intron-3': 7,
            'exon-4': 8,
            'intron-4': 9,
            'exon-5': 10,
            'intron-5': 11,
            'exon-6': 12,
            'intron-6': 13,
            'exon-7': 14,
            'intron-7': 15,
            'exon-8': 16,
            'intron-8': 17,
            'exon-9': 18,
            'intron-9': 19,
            'three_prime_UTR': 20
        },
        'HLA-DQA1':
        {
            'five_prime_UTR': 1,
            'exon-1': 2,
            'intron-1': 3,
            'exon-2': 4,
            'intron-2': 5,
            'exon-3': 6,
            'intron-3': 7,
            'exon-4': 8,
            'three_prime_UTR': 9
        },
        'HLA-DRB3':
        {
            'five_prime_UTR': 1,
            'exon-1': 2,
            'intron-1': 3,
            'exon-2': 4,
            'intron-2': 5,
            'exon-3': 6,
            'intron-3': 7,
            'exon-4': 8,
            'intron-4': 9,
            'exon-5': 10,
            'intron-5': 11,
            'exon-6': 12,
            'three_prime_UTR': 13
        },
        'KIR2DS3':
        {
            'five_prime_UTR': 1,
            'exon-1': 2,
            'intron-1': 3,
            'exon-2': 4,
            'intron-2': 5,
            'exon-3': 6,
            'intron-3': 7,
            'exon-4': 8,
            'intron-4': 9,
            'exon-5': 10,
            'intron-5': 11,
            'exon-6': 12,
            'intron-6': 13,
            'exon-7': 14,
            'intron-7': 15,
            'exon-8': 16,
            'intron-8': 17,
            'exon-9': 18,
            'intron-9': 19,
            'three_prime_UTR': 20
        },
        'KIR3DL1':
        {
            'five_prime_UTR': 1,
            'exon-1': 2,
            'intron-1': 3,
            'exon-2': 4,
            'intron-2': 5,
            'exon-3': 6,
            'intron-3': 7,
            'exon-4': 8,
            'intron-4': 9,
            'exon-5': 10,
            'intron-5': 11,
            'exon-6': 12,
            'intron-6': 13,
            'exon-7': 14,
            'intron-7': 15,
            'exon-8': 16,
            'intron-8': 17,
            'exon-9': 18,
            'intron-9': 19,
            'three_prime_UTR': 20
        },
        'HLA-A':
        {
            'five_prime_UTR': 1,
            'exon-1': 2,
            'intron-1': 3,
            'exon-2': 4,
            'intron-2': 5,
            'exon-3': 6,
            'intron-3': 7,
            'exon-4': 8,
            'intron-4': 9,
            'exon-5': 10,
            'intron-5': 11,
            'exon-6': 12,
            'intron-6': 13,
            'exon-7': 14,
            'intron-7': 15,
            'exon-8': 16,
            'three_prime_UTR': 17
        },
        'HLA-DRB5':
        {
            'five_prime_UTR': 1,
            'exon-1': 2,
            'intron-1': 3,
            'exon-2': 4,
            'intron-2': 5,
            'exon-3': 6,
            'intron-3': 7,
            'exon-4': 8,
            'intron-4': 9,
            'exon-5': 10,
            'intron-5': 11,
            'exon-6': 12,
            'three_prime_UTR': 13
        },
        'KIR2DL4':
        {
            'five_prime_UTR': 1,
            'exon-1': 2,
            'intron-1': 3,
            'exon-2': 4,
            'intron-2': 5,
            'exon-3': 6,
            'intron-3': 7,
            'exon-4': 8,
            'intron-4': 9,
            'exon-5': 10,
            'intron-5': 11,
            'exon-6': 12,
            'intron-6': 13,
            'exon-7': 14,
            'intron-7': 15,
            'exon-8': 16,
            'intron-8': 17,
            'exon-9': 18,
            'intron-9': 19,
            'three_prime_UTR': 20
        },
        'HLA-DQB1':
        {
            'five_prime_UTR': 1,
            'exon-1': 2,
            'intron-1': 3,
            'exon-2': 4,
            'intron-2': 5,
            'exon-3': 6,
            'intron-3': 7,
            'exon-4': 8,
            'intron-4': 9,
            'exon-5': 10,
            'intron-5': 11,
            'exon-6': 12,
            'three_prime_UTR': 13
        },
        'KIR3DL2':
        {
            'five_prime_UTR': 1,
            'exon-1': 2,
            'intron-1': 3,
            'exon-2': 4,
            'intron-2': 5,
            'exon-3': 6,
            'intron-3': 7,
            'exon-4': 8,
            'intron-4': 9,
            'exon-5': 10,
            'intron-5': 11,
            'exon-6': 12,
            'intron-6': 13,
            'exon-7': 14,
            'intron-7': 15,
            'exon-8': 16,
            'intron-8': 17,
            'exon-9': 18,
            'intron-9': 19,
            'three_prime_UTR': 20
        },
        'HLA-B':
        {
            'five_prime_UTR': 1,
            'exon-1': 2,
            'intron-1': 3,
            'exon-2': 4,
            'intron-2': 5,
            'exon-3': 6,
            'intron-3': 7,
            'exon-4': 8,
            'intron-4': 9,
            'exon-5': 10,
            'intron-5': 11,
            'exon-6': 12,
            'intron-6': 13,
            'exon-7': 14,
            'three_prime_UTR': 15
        },
        'KIR3DS1':
        {
            'five_prime_UTR': 1,
            'exon-1': 2,
            'intron-1': 3,
            'exon-2': 4,
            'intron-2': 5,
            'exon-3': 6,
            'intron-3': 7,
            'exon-4': 8,
            'intron-4': 9,
            'exon-5': 10,
            'intron-5': 11,
            'exon-6': 12,
            'intron-6': 13,
            'exon-7': 14,
            'intron-7': 15,
            'exon-8': 16,
            'intron-8': 17,
            'exon-9': 18,
            'intron-9': 19,
            'three_prime_UTR': 20
        },
        'KIR2DL5B':
        {
            'five_prime_UTR': 1,
            'exon-1': 2,
            'intron-1': 3,
            'exon-2': 4,
            'intron-2': 5,
            'exon-3': 6,
            'intron-3': 7,
            'exon-4': 8,
            'intron-4': 9,
            'exon-5': 10,
            'intron-5': 11,
            'exon-6': 12,
            'intron-6': 13,
            'exon-7': 14,
            'intron-7': 15,
            'exon-8': 16,
            'intron-8': 17,
            'exon-9': 18,
            'intron-9': 19,
            'three_prime_UTR': 20
        },
        'HLA-DRB1':
        {
            'five_prime_UTR': 1,
            'exon-1': 2,
            'intron-1': 3,
            'exon-2': 4,
            'intron-2': 5,
            'exon-3': 6,
            'intron-3': 7,
            'exon-4': 8,
            'intron-4': 9,
            'exon-5': 10,
            'intron-5': 11,
            'exon-6': 12,
            'three_prime_UTR': 13
        },
        'KIR3DL3':
        {
            'five_prime_UTR': 1,
            'exon-1': 2,
            'intron-1': 3,
            'exon-2': 4,
            'intron-2': 5,
            'exon-3': 6,
            'intron-3': 7,
            'exon-4': 8,
            'intron-4': 9,
            'exon-5': 10,
            'intron-5': 11,
            'exon-6': 12,
            'intron-6': 13,
            'exon-7': 14,
            'intron-7': 15,
            'exon-8': 16,
            'intron-8': 17,
            'exon-9': 18,
            'intron-9': 19,
            'three_prime_UTR': 20
        },
        'KIR2DS1':
        {
            'five_prime_UTR': 1,
            'exon-1': 2,
            'intron-1': 3,
            'exon-2': 4,
            'intron-2': 5,
            'exon-3': 6,
            'intron-3': 7,
            'exon-4': 8,
            'intron-4': 9,
            'exon-5': 10,
            'intron-5': 11,
            'exon-6': 12,
            'intron-6': 13,
            'exon-7': 14,
            'intron-7': 15,
            'exon-8': 16,
            'intron-8': 17,
            'exon-9': 18,
            'intron-9': 19,
            'three_prime_UTR': 20
        },
        'HLA-DPA1':
        {
            'five_prime_UTR': 1,
            'exon-1': 2,
            'intron-1': 3,
            'exon-2': 4,
            'intron-2': 5,
            'exon-3': 6,
            'intron-3': 7,
            'exon-4': 8,
            'three_prime_UTR': 9
        },
        'HLA-C':
        {
            'five_prime_UTR': 1,
            'exon-1': 2,
            'intron-1': 3,
            'exon-2': 4,
            'intron-2': 5,
            'exon-3': 6,
            'intron-3': 7,
            'exon-4': 8,
            'intron-4': 9,
            'exon-5': 10,
            'intron-5': 11,
            'exon-6': 12,
            'intron-6': 13,
            'exon-7': 14,
            'intron-7': 15,
            'exon-8': 16,
            'three_prime_UTR': 17
        }
    }


