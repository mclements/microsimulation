cat x_gen.c | sed 's/[<]methods\/\(.*\)[>]/"\1"/g' | sed 's/[<]distr\/\(.*\)[>]/"\1"/g' | sed 's/[<]unur[_]source[.]h[>]/"unur_source.h"/g' | head -n 50
