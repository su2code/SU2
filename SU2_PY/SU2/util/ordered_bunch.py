#!/usr/bin/env python

""" OrderedBunch is a subclass of OrderedDict with attribute-style access.

    >>> b = OrderedBunch()
    >>> b.hello = 'world'
    >>> b.hello
    'world'
    >>> b['hello'] += "!"
    >>> b.hello
    'world!'
    >>> b.foo = OrderedBunch(lol=True)
    >>> b.foo.lol
    True
    >>> b.foo is b['foo']
    True

    It is safe to import * from this module:

        __all__ = ('OrderedBunch', 'ordered_bunchify','ordered_unbunchify')

    ordered_un/bunchify provide dictionary conversion; Bunches can also be
    converted via OrderedBunch.to/fromOrderedDict().

    original source:
    https://pypi.python.org/pypi/bunch
"""

from .ordered_dict import OrderedDict

## Compatability Issues...
# try:
#    from collections import OrderedDict
# except ImportError:
#    from ordered_dict import OrderedDict


class OrderedBunch(OrderedDict):
    """A dictionary that provides attribute-style access.

    >>> b = OrderedBunch()
    >>> b.hello = 'world'
    >>> b.hello
    'world'
    >>> b['hello'] += "!"
    >>> b.hello
    'world!'
    >>> b.foo = OrderedBunch(lol=True)
    >>> b.foo.lol
    True
    >>> b.foo is b['foo']
    True

    A OrderedBunch is a subclass of dict; it supports all the methods a dict does...

    >>> b.keys()
    ['foo', 'hello']

    Including update()...

    >>> b.update({ 'ponies': 'are pretty!' }, hello=42)
    >>> print(repr(b))
    OrderedBunch(foo=OrderedBunch(lol=True), hello=42, ponies='are pretty!')

    As well as iteration...

    >>> [ (k,b[k]) for k in b ]
    [('ponies', 'are pretty!'), ('foo', OrderedBunch(lol=True)), ('hello', 42)]

    And "splats".

    >>> "The {knights} who say {ni}!".format(**OrderedBunch(knights='lolcats', ni='can haz'))
    'The lolcats who say can haz!'

    See ordered_unbunchify/OrderedBunch.toOrderedDict, ordered_bunchify/OrderedBunch.fromOrderedDict for notes about conversion.
    """

    _initialized = False

    def __init__(self, *args, **kwarg):
        """initializes the ordered dict"""
        super(OrderedBunch, self).__init__(*args, **kwarg)
        self._initialized = True

    def __contains__(self, k):
        """>>> b = OrderedBunch(ponies='are pretty!')
        >>> 'ponies' in b
        True
        >>> 'foo' in b
        False
        >>> b['foo'] = 42
        >>> 'foo' in b
        True
        >>> b.hello = 'hai'
        >>> 'hello' in b
        True
        """
        try:
            return hasattr(self, k) or dict.__contains__(self, k)
        except:
            return False

    # only called if k not found in normal places
    def __getattr__(self, k):
        """Gets key if it exists, otherwise throws AttributeError.

        nb. __getattr__ is only called if key is not found in normal places.

        >>> b = OrderedBunch(bar='baz', lol={})
        >>> b.foo
        Traceback (most recent call last):
            ...
        AttributeError: foo

        >>> b.bar
        'baz'
        >>> getattr(b, 'bar')
        'baz'
        >>> b['bar']
        'baz'

        >>> b.lol is b['lol']
        True
        >>> b.lol is getattr(b, 'lol')
        True
        """
        try:
            # Throws exception if not in prototype chain
            return object.__getattribute__(self, k)
        except AttributeError:
            try:
                return self[k]
            except KeyError:
                raise AttributeError(k)

    def __setattr__(self, k, v):
        """Sets attribute k if it exists, otherwise sets key k. A KeyError
        raised by set-item (only likely if you subclass OrderedBunch) will
        propagate as an AttributeError instead.

        >>> b = OrderedBunch(foo='bar', this_is='useful when subclassing')
        >>> b.values                            #doctest: +ELLIPSIS
        <built-in method values of OrderedBunch object at 0x...>
        >>> b.values = 'uh oh'
        >>> b.values
        'uh oh'
        >>> b['values']
        Traceback (most recent call last):
            ...
        KeyError: 'values'
        """

        if not self._initialized:
            # for OrderedDict initialization
            return object.__setattr__(self, k, v)

        try:
            # Throws exception if not in prototype chain
            object.__getattribute__(self, k)
        except AttributeError:
            try:
                self[k] = v
            except:
                raise AttributeError(k)
        else:
            object.__setattr__(self, k, v)

    def __delattr__(self, k):
        """Deletes attribute k if it exists, otherwise deletes key k. A KeyError
        raised by deleting the key--such as when the key is missing--will
        propagate as an AttributeError instead.

        >>> b = OrderedBunch(lol=42)
        >>> del b.values
        Traceback (most recent call last):
            ...
        AttributeError: 'OrderedBunch' object attribute 'values' is read-only
        >>> del b.lol
        >>> b.lol
        Traceback (most recent call last):
            ...
        AttributeError: lol
        """
        try:
            # Throws exception if not in prototype chain
            object.__getattribute__(self, k)
        except AttributeError:
            try:
                del self[k]
            except KeyError:
                raise AttributeError(k)
        else:
            object.__delattr__(self, k)

    def toOrderedDict(self):
        """Recursively converts a bunch back into a dictionary.

        >>> b = OrderedBunch(OrderedBunchunch(lol=True), hello=42, ponies='are pretty!')
        >>> b.toOrderedDict()
        {'ponies': 'are pretty!', 'foo': {'lol': True}, 'hello': 42}

        See ordered_unbunchify for more info.
        """
        return ordered_unbunchify(self)

    def __repr__(self):
        """Invertible* string-form of a OrderedBunch.

        >>> b = OrderedBunch(foo=OrderedBunch(lol=True), hello=42, ponies='are pretty!')
        >>> print(repr(b))
        OrderedBunch(foo=OrderedBunch(lol=True), hello=42, ponies='are pretty!')
        >>> eval(repr(b))
        OrderedBunch(foo=OrderedBunch(lol=True), hello=42, ponies='are pretty!')

        (*) Invertible so long as collection contents are each repr-invertible.
        """
        keys = self.keys()
        args = ", ".join(["%s=%r" % (key, self[key]) for key in keys])
        return "%s(%s)" % (self.__class__.__name__, args)

    def __str__(self):
        """String-form of a Bunch."""
        keys = self.keys()
        args = ", ".join(["%s=%r" % (key, self[key]) for key in keys])
        return "{%s}" % args

    @staticmethod
    def fromOrderedDict(d):
        """Recursively transforms a dictionary into a OrderedBunch via copy.

        >>> b = OrderedBunch.fromOrderedDict({'urmom': {'sez': {'what': 'what'}}})
        >>> b.urmom.sez.what
        'what'

        See ordered_bunchify for more info.
        """
        return ordered_bunchify(d)


# While we could convert abstract types like Mapping or Iterable, I think
# ordered_bunchify is more likely to "do what you mean" if it is conservative about
# casting (ex: isinstance(str,Iterable) == True ).
#
# Should you disagree, it is not difficult to duplicate this function with
# more aggressive coercion to suit your own purposes.


def ordered_bunchify(x):
    """Recursively transforms a dictionary into a OrderedBunch via copy.

    >>> b = ordered_bunchify({'urmom': {'sez': {'what': 'what'}}})
    >>> b.urmom.sez.what
    'what'

    ordered_bunchify can handle intermediary dicts, lists and tuples (as well as
    their subclasses), but ymmv on custom datatypes.

    >>> b = ordered_bunchify({ 'lol': ('cats', {'hah':'i win again'}),
    ...         'hello': [{'french':'salut', 'german':'hallo'}] })
    >>> b.hello[0].french
    'salut'
    >>> b.lol[1].hah
    'i win again'

    nb. As dicts are not hashable, they cannot be nested in sets/frozensets.
    """
    if isinstance(x, dict):
        return OrderedBunch((k, ordered_bunchify(v)) for k, v in x.iteritems())
    elif isinstance(x, (list, tuple)):
        return type(x)(ordered_bunchify(v) for v in x)
    else:
        return x


def ordered_unbunchify(x):
    """Recursively converts a OrderedBunch into a dictionary.

    >>> b = OrderedBunch(foo=OrderedBunch(lol=True), hello=42, ponies='are pretty!')
    >>> ordered_unbunchify(b)
    {'ponies': 'are pretty!', 'foo': {'lol': True}, 'hello': 42}

    ordered_unbunchify will handle intermediary dicts, lists and tuples (as well as
    their subclasses), but ymmv on custom datatypes.

    >>> b = OrderedBunch(foo=['bar', OrderedBunch(lol=True)], hello=42,
    ...         ponies=('are pretty!', OrderedBunch(lies='are trouble!')))
    >>> ordered_unbunchify(b) #doctest: +NORMALIZE_WHITESPACE
    {'ponies': ('are pretty!', {'lies': 'are trouble!'}),
     'foo': ['bar', {'lol': True}], 'hello': 42}

    nb. As dicts are not hashable, they cannot be nested in sets/frozensets.
    """
    if isinstance(x, OrderedDict):
        return OrderedDict((k, ordered_unbunchify(v)) for k, v in x.iteritems())
    elif isinstance(x, dict):
        return dict((k, ordered_unbunchify(v)) for k, v in x.iteritems())
    elif isinstance(x, (list, tuple)):
        return type(x)(ordered_unbunchify(v) for v in x)
    else:
        return x


### Serialization

try:
    try:
        import json
    except ImportError:
        import simplejson as json

    def toJSON(self, **options):
        """Serializes this OrderedBunch to JSON. Accepts the same keyword options as `json.dumps()`.

        >>> b = OrderedBunch(foo=OrderedBunch(lol=True), hello=42, ponies='are pretty!')
        >>> json.dumps(b)
        '{"ponies": "are pretty!", "foo": {"lol": true}, "hello": 42}'
        >>> b.toJSON()
        '{"ponies": "are pretty!", "foo": {"lol": true}, "hello": 42}'
        """
        return json.dumps(self, **options)

    OrderedBunch.toJSON = toJSON

except ImportError:
    pass


try:
    # Attempt to register ourself with PyYAML as a representer
    import yaml
    from yaml.representer import Representer, SafeRepresenter

    def from_yaml(loader, node):
        """PyYAML support for Bunches using the tag `!bunch` and `!bunch.OrderedBunch`.

        >>> import yaml
        >>> yaml.load('''
        ... Flow style: !bunch.OrderedBunch { Clark: Evans, Brian: Ingerson, Oren: Ben-Kiki }
        ... Block style: !bunch
        ...   Clark : Evans
        ...   Brian : Ingerson
        ...   Oren  : Ben-Kiki
        ... ''') #doctest: +NORMALIZE_WHITESPACE
        {'Flow style': OrderedBunch(Brian='Ingerson', Clark='Evans', Oren='Ben-Kiki'),
         'Block style': OrderedBunch(Brian='Ingerson', Clark='Evans', Oren='Ben-Kiki')}

        This module registers itself automatically to cover both OrderedBunch and any
        subclasses. Should you want to customize the representation of a subclass,
        simply register it with PyYAML yourself.
        """
        data = OrderedBunch()
        yield data
        value = loader.construct_mapping(node)
        data.update(value)

    def to_yaml_safe(dumper, data):
        """Converts OrderedBunch to a normal mapping node, making it appear as a
        dict in the YAML output.

        >>> b = OrderedBunch(foo=['bar', OrderedBunch(lol=True)], hello=42)
        >>> import yaml
        >>> yaml.safe_dump(b, default_flow_style=True)
        '{foo: [bar, {lol: true}], hello: 42}\\n'
        """
        return dumper.represent_dict(data)

    def to_yaml(dumper, data):
        """Converts OrderedBunch to a representation node.

        >>> b = OrderedBunch(foo=['bar', OrderedBunch(lol=True)], hello=42)
        >>> import yaml
        >>> yaml.dump(b, default_flow_style=True)
        '!bunch.OrderedBunch {foo: [bar, !bunch.OrderedBunch {lol: true}], hello: 42}\\n'
        """
        return dumper.represent_mapping("!orderedbunch.OrderedBunch", data)

    yaml.add_constructor("!orderedbunch", from_yaml)
    yaml.add_constructor("!orderedbunch.OrderedBunch", from_yaml)

    SafeRepresenter.add_representer(OrderedBunch, to_yaml_safe)
    SafeRepresenter.add_multi_representer(OrderedBunch, to_yaml_safe)

    Representer.add_representer(OrderedBunch, to_yaml)
    Representer.add_multi_representer(OrderedBunch, to_yaml)

    # Instance methods for YAML conversion
    def toYAML(self, **options):
        """Serializes this OrderedBunch to YAML, using `yaml.safe_dump()` if
        no `Dumper` is provided. See the PyYAML documentation for more info.

        >>> b = OrderedBunch(foo=['bar', OrderedBunch(lol=True)], hello=42)
        >>> import yaml
        >>> yaml.safe_dump(b, default_flow_style=True)
        '{foo: [bar, {lol: true}], hello: 42}\\n'
        >>> b.toYAML(default_flow_style=True)
        '{foo: [bar, {lol: true}], hello: 42}\\n'
        >>> yaml.dump(b, default_flow_style=True)
        '!bunch.OrderedBunch {foo: [bar, !bunch.OrderedBunch {lol: true}], hello: 42}\\n'
        >>> b.toYAML(Dumper=yaml.Dumper, default_flow_style=True)
        '!bunch.OrderedBunch {foo: [bar, !bunch.OrderedBunch {lol: true}], hello: 42}\\n'
        """
        opts = dict(indent=4, default_flow_style=False)
        opts.update(options)
        if "Dumper" not in opts:
            return yaml.safe_dump(self, **opts)
        else:
            return yaml.dump(self, **opts)

    def fromYAML(*args, **kwargs):
        return ordered_bunchify(yaml.load(*args, **kwargs))

    OrderedBunch.toYAML = OrderedBunch.__repr__ = toYAML
    OrderedBunch.fromYAML = staticmethod(fromYAML)

except ImportError:
    pass
