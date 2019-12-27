import numpy as np
import PIL.Image

def align_lineROI(ref_points, ref_size, img_points, img, bg=None, resample=PIL.Image.NEAREST):
    '''
    Align 2 images by 2 corresponding points (similarity transform: translation, rotation, and scale).
    
    Args
    - ref_points: np.ndarray, shape=(2,2)
        2 points in reference image (0-indexed coordinates, where (0, 0) is the upper left corner).
        x-axis gives column coordinates, y-axis gives row coordinates.
        [[x1, x2],
         [y1, y2]]
    - ref_size: 2-tuple of int
        (width, height)
    - img_points: np.ndarray, shape=(2,2)
        2 points in `img` corresponding to `ref_points`
        [[x1', x2'],
         [y1', y2']]
    - img: PIL.Image.Image
        Image to align to reference.
    - bg: depends. default=None
        Background color used to fill parts of aligned image that are not in the input `img`. If None, defaults to black.
        Value type depends on the mode of the input `img`.
        See https://pillow.readthedocs.io/en/stable/handbook/concepts.html#modes.
    - resample: int. default=0 (PIL.Image.NEAREST)
        Resampling method
    
    Returns: PIL.Image.Image. size=ref_size
      Image aligned to reference image with same size.
    '''
    ref_len = np.sqrt(np.sum((ref_points[:,0] - ref_points[:,1])**2))
    img_len = np.sqrt(np.sum((img_points[:,0] - img_points[:,1])**2))
    scale_factor = ref_len/img_len
    ref_angle = np.arctan2((ref_points[1,1] - ref_points[1,0]), (ref_points[0,1] - ref_points[0,0]))
    img_angle = np.arctan2((img_points[1,1] - img_points[1,0]), (img_points[0,1] - img_points[0,0]))
    theta = ref_angle - img_angle # angle in Cartesian coordinates
    if abs(theta) > np.pi/2:
        raise ValueError(f'The desired rotation ({theta} rad) is beyond pi/2. Please flip the image first if necessary.')
    aligned = img.rotate(-theta * 180 / np.pi, expand=True, fillcolor=bg) # need negative theta to be in image coordinates
    converted = False
    if resample != 0 and aligned.mode in (1, 'P'):
        palette = aligned.getpalette()
        aligned = aligned.convert(mode='L')
        converted = True
    aligned = aligned.resize(np.ceil(np.array(aligned.size) * scale_factor).astype(int), resample=resample)
    if converted:
        aligned = aligned.convert(mode='P', palette=palette)
    
    rot_mat = np.array([[np.cos(theta), -np.sin(theta)],
                        [np.sin(theta),  np.cos(theta)]])
    padding_x = img.size[1] * np.sin(theta) if theta > 0 else 0
    padding_y = img.size[0] * -np.sin(theta) if theta < 0 else 0
    aligned_points = np.round((rot_mat @ img_points + np.array([[padding_x], [padding_y]])) * scale_factor).astype(int)
    
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