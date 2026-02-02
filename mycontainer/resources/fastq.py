class Record:
    # Fields:
    # name
    # description
    # sequence
    # qualities
    def __init__(self, name, description, sequence, qualities, encoding):
        self.name = name
        self.description = description
        self.sequence = sequence
        self.qualities = qualities
        self.encoding = encoding

def Reader(stream, quality_encoding="sanger"):
    while True:
        try:
            id_line = next(stream)
            sequence_line = next(stream)
            id_line2 = next(stream)
            qualities_line = next(stream)
        except StopIteration:
            break
        
        assert(id_line.startswith('@'))
        assert(id_line2.startswith('+'))

        id_fields = id_line[1:-1].split(None, 1)
        name = id_fields[0]
        description = id_fields[1] if len(id_fields) == 2 else ""
        sequence = sequence_line[:-1]
        qualities = qualities_line[:-1]

        yield Record(name, description, sequence, qualities, quality_encoding)

class Writer:
    def __init__(self, stream, quality_encoding="sanger"):
        self.stream = stream

    def write(self, rec):
        title = rec.name
        if rec.description:
            title += " " + rec.description
        self.stream.write("@%s\n%s\n+\n%s\n" % (title,
                                                rec.sequence,
                                                rec.qualities))
