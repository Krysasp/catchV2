# Obfuscated CATCHv2 script
# Original: design_large.py
# Functions: ...

#!/usr/bin/env python3
"""Design probes for genome capture, with options and parameters that
optimize resource usage for large, highly diverse input.

The downside of these options is usually a small increase in the number of
designed probes. As with design.py, this program still supports full
customization of the argument values.

This wraps design.py and offers a way to run design.py without requiring
deep familiarity with CATCH'𐁳𐀠𐁯𐁰𐁴𐁩𐁯𐁮𐁳𐀮𐀠𐁔𐁨𐁡𐁴𐀠𐁩𐁳𐀬𐀠𐁩𐁴𐀠𐁴𐁡𐁫𐁥𐁳𐀠𐁩𐁮𐁴𐁯𐀠𐁡𐁣𐁣𐁯𐁵𐁮𐁴𐀊𐁲𐁥𐁣𐁯𐁭𐁭𐁥𐁮𐁤𐁡𐁴𐁩𐁯𐁮𐁳𐀠𐁴𐁨𐁡𐁴𐀠𐁯𐁦𐁴𐁥𐁮𐀠𐁷𐁯𐁲𐁫𐀠𐁷𐁥𐁬𐁬𐀠𐁩𐁮𐀠𐁰𐁲𐁡𐁣𐁴𐁩𐁣𐁥𐀮𐀊𐀢𐀢𐀢𐀊𐀊𐁩𐁭𐁰𐁯𐁲𐁴𐀠𐁤𐁥𐁳𐁩𐁧𐁮𐀊𐀊𐁟𐁟𐁡𐁵𐁴𐁨𐁯𐁲𐁟𐁟𐀠𐀽𐀠'Hayden Metsky <hayden@broadinstitute.org>'𐀊𐀊𐀊𐁩𐁦𐀠𐁟𐁟𐁮𐁡𐁭𐁥𐁟𐁟𐀠𐀽𐀽𐀠𐀢𐁟𐁟𐁭𐁡𐁩𐁮𐁟𐁟𐀢𐀺𐀊𐀠𐀠𐀠𐀠𐁡𐁲𐁧𐁳𐀠𐀽𐀠𐁤𐁥𐁳𐁩𐁧𐁮𐀮𐁩𐁮𐁩𐁴𐁟𐁡𐁮𐁤𐁟𐁰𐁡𐁲𐁳𐁥𐁟𐁡𐁲𐁧𐁳𐀨𐁡𐁲𐁧𐁳𐁟𐁴𐁹𐁰𐁥𐀽'large')
    design.main(args)
