import warnings
import numpy as np
import PIL.Image

def align_auto(ref_points, ref_size, img_points, img, ref_aspect=1, bg=None, resample=PIL.Image.NEAREST, return_transform=True):
    '''
    Align 2 images by given a set of corresponding points.
    - 2 corresponding points provided --> align_lineROI()
    - 3+ corresponding points provided --> align_affine()

    Args
    - ref_points: np.ndarray, shape=(2,2)
        2 points in reference image (0-indexed coordinates, where (0, 0) is the upper left corner).
        x-axis gives column coordinates, y-axis gives row coordinates.
        [[x1, x2],
         [y1, y2]]
    - ref_size: 2-tuple of int
        Size of reference image: (width, height)
    - img_points: np.ndarray, shape=(2,2)
        2 points in `img` corresponding to `ref_points`
        [[x1', x2'],
         [y1', y2']]
    - img: PIL.Image.Image
        Image to align to reference.
    - ref_aspect: float. default=1
        (inches wide per pixel) / (inches tall per pixel) of reference image, relative to `img`
    - bg: depends. default=None
        Background color used to fill parts of aligned image that are not in the input `img`. If None, defaults to black.
        Value type depends on the mode of the input `img`.
        See https://pillow.readthedocs.io/en/stable/handbook/concepts.html#modes.
    - resample: int. default=PIL.Image.NEAREST
        Resampling method: PIL.Image.(NEAREST|BOX|BILINEAR|HAMMING|BICUBIC|LANCZOS)
        See https://pillow.readthedocs.io/en/stable/handbook/concepts.html#filters.
    - return_transform: bool. default=True
        Return transform, if applicable.

    Returns: PIL.Image.Image. size=ref_size
      Image aligned to reference image with same size.
    '''
    n_points = ref_points.shape[1]
    if n_points == 2:
        aligned, transform = align_lineROI(ref_points, ref_size, img_points, img, ref_aspect, bg, resample), None
    else:
        aligned, transform = align_affine(ref_points, ref_size, img_points, img, ref_aspect, bg, resample)
    if return_transform:
        return aligned, transform
    return aligned

def align_lineROI(ref_points, ref_size, img_points, img, ref_aspect=1, bg=None, resample=PIL.Image.NEAREST):
    '''
    Align 2 images by 2 corresponding points (similarity transform: translation, rotation, and scale).

    Args: see align_auto()

    Returns: PIL.Image.Image. size=ref_size
      Image aligned to reference image.
    '''
    assert ref_points.shape == img_points.shape, 'ref_points and img_points should have the same shape.'
    if ref_points.shape[1] > 2:
        warnings.warn(f'{ref_points.shape[1]} corresponding points provided, but align_lineROI() will only use the first two.')
    ref_points = ref_points.copy()
    if ref_aspect != 1:
        ref_points_orig = ref_points.copy()
        ref_points[1, :] = ref_points[1, :] / ref_aspect   # ref_aspect < 1: each pixel represents more height than suggested

    ref_len = np.sqrt(np.sum((ref_points[:,0] - ref_points[:,1])**2))
    img_len = np.sqrt(np.sum((img_points[:,0] - img_points[:,1])**2))
    scale_factor = ref_len/img_len

    # rotate
    # - PIL.Image.rotate(deg) rotates the image `deg` degrees counterclockwise around its center
    ref_angle = np.arctan2((ref_points[1,1] - ref_points[1,0]), (ref_points[0,1] - ref_points[0,0]))
    img_angle = np.arctan2((img_points[1,1] - img_points[1,0]), (img_points[0,1] - img_points[0,0]))
    theta = ref_angle - img_angle # radians; positive value: rotate img clockwise to match orientation of ref
    if abs(theta) > np.pi/2:
        raise ValueError((f'The magnitude of the desired rotation ({theta} rad) is beyond pi/2. '
                          'Please flip the image first if necessary.'))
    aligned = img.rotate(-theta * 180 / np.pi, expand=True, fillcolor=bg)

    # resize
    # - if img is black-and-white (mode=='1') or an 8-bit image with custom color palette,
    #   and we want to use an interpolating resampling method during resizing,
    #   convert img to grayscale before resizing
    converted = False
    if resample != PIL.Image.NEAREST and aligned.mode in ('1', 'P'):
        palette = aligned.getpalette()
        aligned = aligned.convert(mode='L')
        converted = True
    width_new = int(round(aligned.size[0] * scale_factor))
    height_new = int(round(aligned.size[1] * scale_factor * ref_aspect))
    aligned = aligned.resize((width_new, height_new), resample=resample)
    if converted:
        aligned = aligned.convert(mode='P', palette=palette)

    rot_mat = np.array([[np.cos(theta), -np.sin(theta)],
                        [np.sin(theta),  np.cos(theta)]])
    padding_x = img.size[1] * np.sin(theta) if theta > 0 else 0
    padding_y = img.size[0] * -np.sin(theta) if theta < 0 else 0
    aligned_points = (rot_mat @ img_points + np.array([[padding_x], [padding_y]])) * scale_factor
    if ref_aspect != 1:
        aligned_points[1, :] *= ref_aspect
        ref_points = ref_points_orig.copy()
    aligned_points = np.round(aligned_points).astype(int)

    crop_box = (
        aligned_points[0,0] - ref_points[0,0],               # left
        aligned_points[1,0] - ref_points[1,0],               # upper
        aligned_points[0,0] - ref_points[0,0] + ref_size[0], # right
        aligned_points[1,0] - ref_points[1,0] + ref_size[1]  # lower
    )

    if (any(np.array(crop_box) < 0) or crop_box[2] > aligned.size[0] or crop_box[3] > aligned.size[1]) and bg is not None:
        expand_x = abs(min(crop_box[0], 0))
        expand_y = abs(min(crop_box[1], 0))
        tmp = PIL.Image.new(img.mode, (ref_size[0] + expand_x, ref_size[1] + expand_y), color=bg)
        tmp.paste(aligned, (expand_x, expand_y))
        return tmp
    else:
        return aligned.crop(crop_box)

def align_affine(ref_points, ref_size, img_points, img, ref_aspect=1, bg=None, resample=PIL.Image.NEAREST):
    '''
    Aligns 2 images using an affine transformation given 3+ corresponding points.

    Args: see align_auto()

    Returns
    - aligned: PIL.Image.Image. size=ref_size
        Image aligned to reference image.
    - transform: numpy.ndarray. shape=(3, 3)
        Affine transformation matrix
    '''
    n_points = ref_points.shape[1]
    assert n_points >= 3, 'Need at least 3 points to define affine transformation.'
    assert ref_points.shape == img_points.shape, 'ref_points and img_points should have the same shape.'
    if ref_points.shape[0] == 3:
        assert np.all(ref_points[2, :] == 1) and np.all(img_points[2, :] == 1)
    else:
        assert ref_points.shape[0] == 2
        ref_points = np.vstack((ref_points, np.ones(n_points)))
        img_points = np.vstack((img_points, np.ones(n_points)))

    ref_points = ref_points.copy()
    ref_points_orig = ref_points.copy()
    ref_points[1, :] = ref_points[1, :] / ref_aspect   # ref_aspect < 1: each pixel represents more height than suggested
    scale = np.eye(3)
    scale[1, 1] /= ref_aspect
    transform = scale @ img_points @ np.linalg.pinv(ref_points)

    aligned = img.transform(
        ref_size,
        method=PIL.Image.AFFINE,
        data=transform.ravel()[:6],
        fillcolor=bg,
        resample=resample)
    # for images represented as numpy arrays rather than pillow images:
    # skimage.transform.warp(np.array(img), skimage.transform.AffineTransform(transform), output_shape=ref_size, order=3)
    return aligned, transform
