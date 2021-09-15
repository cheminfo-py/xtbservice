from diskcache import Cache

opt_cache = Cache(size_limit=2 * 10 ** 7, disk_min_file_size=0)
ir_cache = Cache(size_limit=2 * 10 ** 7, disk_min_file_size=0)