#!/usr/bin/env python

""" Bunch is a subclass of dict with attribute-style access.

    >>> b = Bunch()
    >>> b.hello = 'world'
    >>> b.hello
    'world'
    >>> b['hello'] += "!"
    >>> b.hello
    'world!'
    >>> b.foo = Bunch(lol=True)
    >>> b.foo.lol
    True
    >>> b.foo is b['foo']
    True

    It is safe to import * from this module:

        __all__ = ('Bunch', 'bunchify','unbunchify')

    un/bunchify provide dictionary conversion; Bunches can also be
    converted via Bunch.to/fromDict().

    original source:
    https://pypi.python.org/pypi/bunch
"""


class Bunch(dict):
    """A dictionary that provides attribute-style access.

    >>> b = Bunch()
    >>> b.hello = 'world'
    >>> b.hello
    'world'
    >>> b['hello'] += "!"
    >>> b.hello
    'world!'
    >>> b.foo = Bunch(lol=True)
    >>> b.foo.lol
    True
    >>> b.foo is b['foo']
    True

    A Bunch is a subclass of dict; it supports all the methods a dict does...

    >>> b.keys()
    ['foo', 'hello']

    Including update()...

    >>> b.update({ 'ponies': 'are pretty!' }, hello=42)
    >>> print(repr(b))
    Bunch(foo=Bunch(lol=True), hello=42, ponies='are pretty!')

    As well as iteration...

    >>> [ (k,b[k]) for k in b ]
    [('ponies', 'are pretty!'), ('foo', Bunch(lol=True)), ('hello', 42)]

    And "splats".

    >>> "The {knights} who say {ni}!".format(**Bunch(knights='lolcats', ni='can haz'))
    'The lolcats who say can haz!'

    See unbunchify/Bunch.toDict, bunchify/Bunch.fromDict for notes about conversion.
    """

    def __contains__(self, k):
        """>>> b = Bunch(ponies='are pretty!')
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

        >>> b = Bunch(bar='baz', lol={})
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
        raised by set-item (only likely if you subclass Bunch) will
        propagate as an AttributeError instead.

        >>> b = Bunch(foo='bar', this_is='useful when subclassing')
        >>> b.values                            #doctest: +ELLIPSIS
        <built-in method values of Bunch object at 0x...>
        >>> b.values = 'uh oh'
        >>> b.values
        'uh oh'
        >>> b['values']
        Traceback (most recent call last):
            ...
        KeyError: 'values'
        """
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

        >>> b = Bunch(lol=42)
        >>> del b.values
        Traceback (most recent call last):
            ...
        AttributeError: 'Bunch' object attribute 'values' is read-only
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

    def toDict(self):
        """Recursively converts a bunch back into a dictionary.

        >>> b = Bunch(foo=Bunch(lol=True), hello=42, ponies='are pretty!')
        >>> b.toDict()
        {'ponies': 'are pretty!', 'foo': {'lol': True}, 'hello': 42}

        See unbunchify for more info.
        """
        return unbunchify(self)

    def __repr__(self):
        """Invertible* string-form of a Bunch.

        >>> b = Bunch(foo=Bunch(lol=True), hello=42, ponies='are pretty!')
        >>> print(repr(b))
        Bunch(foo=Bunch(lol=True), hello=42, ponies='are pretty!')
        >>> eval(repr(b))
        Bunch(foo=Bunch(lol=True), hello=42, ponies='are pretty!')

        (*) Invertible so long as collection contents are each repr-invertible.
        """
        keys = self.keys()
        keys.sort()
        args = ", ".join(["%s=%r" % (key, self[key]) for key in keys])
        return "%s(%s)" % (self.__class__.__name__, args)

    def __str__(self):
        """String-form of a OrderedBunch."""
        keys = self.keys()
        keys.sort()
        args = ", ".join(["%s=%r" % (key, self[key]) for key in keys])
        return "{%s}" % args

    @staticmethod
    def fromDict(d):
        """Recursively transforms a dictionary into a Bunch via copy.

        >>> b = Bunch.fromDict({'urmom': {'sez': {'what': 'what'}}})
        >>> b.urmom.sez.what
        'what'

        See bunchify for more info.
        """
        return bunchify(d)


# While we could convert abstract types like Mapping or Iterable, I think
# bunchify is more likely to "do what you mean" if it is conservative about
# casting (ex: isinstance(str,Iterable) == True ).
#
# Should you disagree, it is not difficult to duplicate this function with
# more aggressive coercion to suit your own purposes.


def bunchify(x):
    """Recursively transforms a dictionary into a Bunch via copy.

    >>> b = bunchify({'urmom': {'sez': {'what': 'what'}}})
    >>> b.urmom.sez.what
    'what'

    bunchify can handle intermediary dicts, lists and tuples (as well as
    their subclasses), but ymmv on custom datatypes.

    >>> b = bunchify({ 'lol': ('cats', {'hah':'i win again'}),
    ...         'hello': [{'french':'salut', 'german':'hallo'}] })
    >>> b.hello[0].french
    'salut'
    >>> b.lol[1].hah
    'i win again'

    nb. As dicts are not hashable, they cannot be nested in sets/frozensets.
    """
    if isinstance(x, dict):
        return Bunch((k, bunchify(v)) for k, v in x.iteritems())
    elif isinstance(x, (list, tuple)):
        return type(x)(bunchify(v) for v in x)
    else:
        return x


def unbunchify(x):
    """Recursively converts a Bunch into a dictionary.

    >>> b = Bunch(foo=Bunch(lol=True), hello=42, ponies='are pretty!')
    >>> unbunchify(b)
    {'ponies': 'are pretty!', 'foo': {'lol': True}, 'hello': 42}

    unbunchify will handle intermediary dicts, lists and tuples (as well as
    their subclasses), but ymmv on custom datatypes.

    >>> b = Bunch(foo=['bar', Bunch(lol=True)], hello=42,
    ...         ponies=('are pretty!', Bunch(lies='are trouble!')))
    >>> unbunchify(b) #doctest: +NORMALIZE_WHITESPACE
    {'ponies': ('are pretty!', {'lies': 'are trouble!'}),
     'foo': ['bar', {'lol': True}], 'hello': 42}

    nb. As dicts are not hashable, they cannot be nested in sets/frozensets.
    """
    if isinstance(x, dict):
        return dict((k, unbunchify(v)) for k, v in x.iteritems())
    elif isinstance(x, (list, tuple)):
        return type(x)(unbunchify(v) for v in x)
    else:
        return x


### Serialization

try:
    try:
        import json
    except ImportError:
        import simplejson as json

    def toJSON(self, **options):
        """Serializes this Bunch to JSON. Accepts the same keyword options as `json.dumps()`.

        >>> b = Bunch(foo=Bunch(lol=True), hello=42, ponies='are pretty!')
        >>> json.dumps(b)
        '{"ponies": "are pretty!", "foo": {"lol": true}, "hello": 42}'
        >>> b.toJSON()
        '{"ponies": "are pretty!", "foo": {"lol": true}, "hello": 42}'
        """
        return json.dumps(self, **options)

    Bunch.toJSON = toJSON

except ImportError:
    pass


try:
    # Attempt to register ourself with PyYAML as a representer
    import yaml
    from yaml.representer import Representer, SafeRepresenter

    def from_yaml(loader, node):
        """PyYAML support for Bunches using the tag `!bunch` and `!bunch.Bunch`.

        >>> import yaml
        >>> yaml.load('''
        ... Flow style: !bunch.Bunch { Clark: Evans, Brian: Ingerson, Oren: Ben-Kiki }
        ... Block style: !bunch
        ...   Clark : Evans
        ...   Brian : Ingerson
        ...   Oren  : Ben-Kiki
        ... ''') #doctest: +NORMALIZE_WHITESPACE
        {'Flow style': Bunch(Brian='Ingerson', Clark='Evans', Oren='Ben-Kiki'),
         'Block style': Bunch(Brian='Ingerson', Clark='Evans', Oren='Ben-Kiki')}

        This module registers itself automatically to cover both Bunch and any
        subclasses. Should you want to customize the representation of a subclass,
        simply register it with PyYAML yourself.
        """
        data = Bunch()
        yield data
        value = loader.construct_mapping(node)
        data.update(value)

    def to_yaml_safe(dumper, data):
        """Converts Bunch to a normal mapping node, making it appear as a
        dict in the YAML output.

        >>> b = Bunch(foo=['bar', Bunch(lol=True)], hello=42)
        >>> import yaml
        >>> yaml.safe_dump(b, default_flow_style=True)
        '{foo: [bar, {lol: true}], hello: 42}\\n'
        """
        return dumper.represent_dict(data)

    def to_yaml(dumper, data):
        """Converts Bunch to a representation node.

        >>> b = Bunch(foo=['bar', Bunch(lol=True)], hello=42)
        >>> import yaml
        >>> yaml.dump(b, default_flow_style=True)
        '!bunch.Bunch {foo: [bar, !bunch.Bunch {lol: true}], hello: 42}\\n'
        """
        return dumper.represent_mapping("!bunch.Bunch", data)

    yaml.add_constructor("!bunch", from_yaml)
    yaml.add_constructor("!bunch.Bunch", from_yaml)

    SafeRepresenter.add_representer(Bunch, to_yaml_safe)
    SafeRepresenter.add_multi_representer(Bunch, to_yaml_safe)

    Representer.add_representer(Bunch, to_yaml)
    Representer.add_multi_representer(Bunch, to_yaml)

    # Instance methods for YAML conversion
    def toYAML(self, **options):
        """Serializes this Bunch to YAML, using `yaml.safe_dump()` if
        no `Dumper` is provided. See the PyYAML documentation for more info.

        >>> b = Bunch(foo=['bar', Bunch(lol=True)], hello=42)
        >>> import yaml
        >>> yaml.safe_dump(b, default_flow_style=True)
        '{foo: [bar, {lol: true}], hello: 42}\\n'
        >>> b.toYAML(default_flow_style=True)
        '{foo: [bar, {lol: true}], hello: 42}\\n'
        >>> yaml.dump(b, default_flow_style=True)
        '!bunch.Bunch {foo: [bar, !bunch.Bunch {lol: true}], hello: 42}\\n'
        >>> b.toYAML(Dumper=yaml.Dumper, default_flow_style=True)
        '!bunch.Bunch {foo: [bar, !bunch.Bunch {lol: true}], hello: 42}\\n'
        """
        opts = dict(indent=4, default_flow_style=False)
        opts.update(options)
        if "Dumper" not in opts:
            return yaml.safe_dump(self, **opts)
        else:
            return yaml.dump(self, **opts)

    def fromYAML(*args, **kwargs):
        return bunchify(yaml.load(*args, **kwargs))

    Bunch.toYAML = Bunch.__repr__ = toYAML
    Bunch.fromYAML = staticmethod(fromYAML)

except ImportError:
    pass
