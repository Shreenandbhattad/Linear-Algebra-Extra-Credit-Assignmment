# Recognition, Detection & Neural Radiance Field Reconstruction

> **Computer Vision Final Project** — Classification · Detection · BEV Transform · NeRF  
> Dataset: [Khana](https://khana.omkar.xyz/) · 131,000+ images · 80 Indian food classes

---

## Table of Contents
1. [Project Overview](#project-overview)
2. [Problem 1 — Image Classification](#problem-1--image-classification)
3. [Problem 2 — Thali Food Detection](#problem-2--thali-food-detection)
4. [Problem 3 — Bird's Eye View Transform](#problem-3--birds-eye-view-transform)
5. [Bonus — Neural Radiance Field](#bonus--neural-radiance-field)
6. [Results Summary](#results-summary)
7. [Repository Structure](#repository-structure)
8. [Setup & Run](#setup--run)

---

## Project Overview

This project addresses four computer vision challenges using the **Khana dataset** — a comprehensive Indian cuisine benchmark containing 131,000+ images across 80 food classes. The dataset was specifically designed to address the lack of representation of Indian food in existing benchmarks such as Food-101, which is heavily Western-biased.

The four tasks span the core pillars of modern computer vision: **image classification**, **instance segmentation-based object detection**, **geometric image transformation**, and **novel view synthesis via implicit neural representations**.

---

## Problem 1 — Image Classification

### Task
Train a classifier on the Khana dataset to exceed the baseline validation accuracy of **~91%** on an 80/20 train-validation split across 80 food classes.

### Architecture
We use **ConvNeXt-Base** (`convnext_base.fb_in22k_ft_in1k`) pretrained on ImageNet-21K and fine-tuned on ImageNet-1K via the `timm` library. ConvNeXt modernises the standard ResNet design by adopting transformer-inspired architectural choices (depthwise convolutions, inverted bottlenecks, GELU activations, LayerNorm) while retaining the inductive biases of CNNs.

### Training Strategy — Two-Stage Progressive Resolution

**Stage 1 (224px, 20 epochs)**
- Optimizer: AdamW, base LR = 3×10⁻⁴, weight decay = 0.05
- Scheduler: Cosine annealing with 5-epoch linear warmup
- Layer-wise LR decay (factor 0.75 per stage) — head gets full LR, stem gets lowest LR
- Stochastic depth (drop path rate = 0.2) for regularisation

**Stage 2 (320px, 8 epochs)**  
- Load best Stage 1 checkpoint, fine-tune at higher resolution
- LR = 5×10⁻⁵ with 1-epoch warmup
- Same cosine schedule

### Data Augmentation
- **RandAugment** (`rand-m9-mstd0.5-inc1`) — magnitude 9, 14 augmentation ops
- **Random Erasing** (probability 0.25, pixel-fill mode)
- **MixUp** (α = 0.4) + **CutMix** (α = 1.0), switching probability 0.5
- Random horizontal flip, colour jitter (strength 0.4), scale jitter (0.65–1.0)

### Class Imbalance Handling
- **WeightedRandomSampler** — each class sampled equally per epoch regardless of image count
- **Soft-target Cross-Entropy** loss (compatible with MixUp/CutMix soft labels)
- Label smoothing = 0.1

### Evaluation
- **Test-Time Augmentation (TTA)**: 5-crop (4 corners + centre), logits averaged
- Metric: Top-1 and Top-5 validation accuracy

### Results

| Metric | Value |
|--------|-------|
| Baseline (given) | ~91.00% |
| Our Val Acc @1 | **95.87%** |
| Our Val Acc @5 | >99% |
| Above baseline | **+4.87%** |
| Architecture | ConvNeXt-Base (88M params) |
| Training resolution | 224px → 320px |

---

## Problem 2 — Thali Food Detection

### Task
Given a clean bird's-eye-view thali image, detect individual food portions and label them with the correct Khana class names — replacing the generic COCO labels (bowl, dining table) produced by off-the-shelf detectors. Evaluated by **Precision/Recall** on labels only (bounding box IoU not evaluated).

### Why Not a Standard Object Detector?
Standard detectors (YOLO, Faster-RCNN) trained on COCO output labels like `bowl`, `dining table`, `sandwich` — not `chana masala`, `chapati`, `biryani`. Fine-tuning a detector from scratch on Khana would require bounding box annotations, which the dataset does not provide (classification only). We instead build a **segmentation-then-classify** pipeline.

### Pipeline

```
Input Image
    │
    ▼
SAM (Segment Anything Model) — vit_h
    │  Automatic mask generation
    │  Points per side: 48, IoU thresh: 0.88, Stability thresh: 0.95
    │
    ▼
Mask Filtering
    │  Area fraction: 1% – 18% of image
    │  Aspect ratio: 0.20 – 5.0
    │  Min crop size: 24×24 px
    │
    ▼
ConvNeXt-Base Classifier (our trained model)
    │  Each mask region cropped + classified
    │  Confidence threshold: 0.75
    │
    ▼
Non-Maximum Suppression (NMS)
    │  Per-class NMS then global NMS
    │  IoU threshold: 0.25
    │
    ▼
Output: Bounding boxes + Khana labels + confidence scores
```

### Key Components

**Segment Anything Model (SAM)**  
SAM (`vit_h`, 2.39 GB) uses a Vision Transformer image encoder, a prompt encoder, and a mask decoder. In automatic mode it densely samples point prompts across the image and generates all plausible masks without any class supervision. This gives us class-agnostic food region proposals.

**Crop Classification**  
Each SAM mask's bounding box is padded by 10px, cropped from the original image, and passed through our trained ConvNeXt-Base classifier. The softmax probability of the top-1 prediction becomes the detection confidence. Top-5 predictions are also stored.

**NMS**  
Custom two-pass NMS: first per-class (removes duplicates of the same food), then global (removes cross-class overlapping boxes). IoU computed as intersection over union of axis-aligned bounding boxes.

### Results (IMG_0901.jpg, BEV mode)

| # | Label | Confidence |
|---|-------|-----------|
| 1 | chapati | 0.912 |
| 2 | chana masala | 0.901 |
| 3 | chana masala | 0.898 |
| 4 | biryani (rice) | 0.755 |

---

## Problem 3 — Bird's Eye View Transform

### Task
Natural thali photos are taken from oblique angles, causing perspective distortion that degrades SAM segmentation quality. Transform natural images to a top-down Bird's Eye View (BEV) before running the detection pipeline.

### Methods Implemented

**Interactive BEV (Primary)**  
User clicks 4 corners of the thali tray in order (TL → TR → BR → BL). A perspective homography `H ∈ ℝ³ˣ³` is computed via `cv2.getPerspectiveTransform` and applied with `cv2.warpPerspective` using bicubic interpolation. Output is square with side = min(H, W) of original.

**Automatic BEV (Fallback)**  
Hough Circle Transform (`cv2.HoughCircles`) detects the circular plate boundary. The bounding square of the detected circle is used as the four source points for the homography. Works well for round plates; less reliable for rectangular thali trays.

### Homography Estimation
Given source points **p** = {p₁, p₂, p₃, p₄} and destination points **p'** forming a square, the homography satisfies:

```
p'ᵢ = H · pᵢ   (in homogeneous coordinates)
```

`cv2.getPerspectiveTransform` solves the 8 degrees of freedom (up to scale) using the 4 point correspondences exactly via direct linear transform (DLT).

### Pipeline
```
Natural Image (oblique angle)
    │
    ▼
BEV Warp (homography)
    │
    ▼
SAM Segmentation
    │
    ▼
ConvNeXt Classification
    │
    ▼
NMS + Output
```

### Qualitative Results
After BEV warping, SAM detects significantly more masks (83 vs 37 on the same image) because food regions become more regular and top-facing, improving the stability score of mask candidates. Detection labels are consistently more accurate post-warp.

---

## Bonus — Neural Radiance Field

### Task
Reconstruct a Neural Radiance Field (NeRF) of a physical thali following the CS180 Project 4 pipeline. This requires camera calibration, multi-view pose estimation, a 2D neural field warmup, and full 3D NeRF training.

---

### Step 0a — Camera Calibration

**Method**: Multi-view ArUco marker calibration using a printed 2×3 grid (IDs 0–5, DICT_4X4_50).

**Board specifications** (from printed PDF):
- Marker size: 60mm
- Horizontal spacing: 90.00mm centre-to-centre
- Vertical spacing: 75.67mm centre-to-centre
- Grid layout: 2 columns × 3 rows

**3D world coordinates** are computed per tag using exact physical spacings (non-square grid — this was the key fix that brought RMS from 72px down to 1.2px). Each tag corner's 3D position:

```
cx = col × 0.090,   cy = row × 0.07567
corners = [cx±s, cy±s, 0]   where s = 0.030
```

`cv2.calibrateCamera` solves for the intrinsic matrix **K** and distortion coefficients **d** using all detected corner correspondences across 57 images.

**Calibration Results**:

| Parameter | Value |
|-----------|-------|
| RMS reprojection error | **1.21 px** |
| fx | 680.5 px |
| fy | 684.5 px |
| cx | 363.6 px |
| cy | 310.2 px |
| fx/fy ratio | 0.994 (expected ~1.0) ✓ |

---

### Step 0b — Camera Pose Estimation

For each thali image, `cv2.solvePnP` (IPPE_SQUARE method) estimates the camera extrinsics relative to the ArUco marker using the calibrated **K** matrix.

**World-to-Camera → Camera-to-World conversion**:
```python
R, _ = cv2.Rodrigues(rvec)
c2w[:3,:3] = R.T
c2w[:3, 3] = -(R.T @ tvec)
```

**Coordinate convention flip** (OpenCV → OpenGL):
```python
c2w[:3, 1] *= -1   # flip Y
c2w[:3, 2] *= -1   # flip Z
```

**Dataset**: 51 of 55 images successfully posed. Split: 40 train / 5 val / 6 test. Images resized to 400×300, focal length: 333.6px.

---

### Step 1 — 2D Neural Field

A 2D MLP is trained to map pixel coordinates (u, v) ∈ [0,1]² to RGB colour. This validates the positional encoding and training setup before full 3D training.

**Architecture**:
- Positional encoding: L=10 frequency levels → input dim = 2 + 4×10 = 42
- MLP: 4 hidden layers × 256 units, ReLU activations, Sigmoid output
- Loss: MSE on RGB values

**Result**: **34.06 dB PSNR** after 3000 iterations (~4 min on CPU)

---

### Step 2 — 3D NeRF

Full volumetric scene representation following Mildenhall et al. (2020).

**Architecture**: Coarse + Fine NeRF MLPs (separate weights)
- Position encoding: L_pos=10 → 63-dim input
- Direction encoding: L_dir=4 → 27-dim input
- 8-layer MLP, hidden dim 256, skip connection at layer 4
- View-dependent colour head: 256+27 → 128 → 3 (RGB)
- Density head: 256 → 1 (σ)

**Rendering**:

Coarse: stratified sampling of N_c=96 points per ray in [near, far]

Fine: hierarchical importance sampling of N_f=192 additional points using coarse weights as PDF

Volume rendering integral (discretised):
```
C(r) = Σᵢ Tᵢ · αᵢ · cᵢ
Tᵢ = exp(-Σⱼ<ᵢ σⱼδⱼ)
αᵢ = 1 - exp(-σᵢδᵢ)
```

**Key improvements over vanilla NeRF**:

1. **Scene normalisation** — all camera positions centred and scaled to unit sphere. NEAR/FAR computed adaptively from pose distribution (NEAR=0.01, FAR=2.41). This fixed the 4.5dB → 22dB improvement.

2. **White background blending** — unoccupied rays rendered as white (matching table background):
   ```
   C_final = C_nerf + (1 - Σ weights)
   ```

3. **Mixed precision training** (AMP) — `torch.amp.autocast('cuda')` + `GradScaler` for ~30% speedup on T4

4. **Gradient clipping** — `clip_grad_norm_(params, 1.0)` for training stability

5. **Cosine LR decay** — LR 5×10⁻⁴ → 5×10⁻⁶ over 20k steps

**Training**: 20,000 iterations, batch 4096 rays, ~2 hours on NVIDIA T4 (Google Colab)

**NeRF Results**:

| Metric | Value |
|--------|-------|
| Final PSNR | **22.00 dB** |
| Training steps | 20,000 |
| Hardware | NVIDIA T4 (16GB) |
| Spiral video | 60 frames @ 24fps |
| Scene scale | 0.47m (normalised) |

---

## Results Summary

| Problem | Method | Metric | Result |
|---------|--------|--------|--------|
| 1. Classification | ConvNeXt-Base + 2-stage training | Val Acc @1 | **95.87%** (+4.87% over baseline) |
| 2. Detection | SAM vit_h + ConvNeXt classifier | Precision/Recall on labels | Correct Khana labels, conf > 0.75 |
| 3. BEV | Perspective homography (4-point DLT) | Qualitative | Clean top-down view, improved SAM masks |
| 4. NeRF (bonus) | Coarse-fine NeRF + scene normalisation | PSNR | **22.00 dB**, spiral GIF produced |

---

## Repository Structure

```
Computer-Vision/
├── codes/
│   ├── train.py                    # Problem 1 - classifier training
│   ├── predict_khana.py            # Problem 1 - inference on test images
│   ├── detect_thali.py             # Problem 2+3 - SAM detection + BEV
│   ├── step0_calibrate.py          # Bonus - camera calibration
│   ├── step0_pose_and_dataset.py   # Bonus - pose estimation + dataset
│   ├── step1_neural_field_2d.py    # Bonus - 2D neural field
│   ├── step2_nerf_3d.py            # Bonus - full 3D NeRF
│   ├── check_results.py            # summary of all results
│   └── test_gpu.py                 # GPU check
├── nerf_results/
│   ├── spiral.gif                  # 60-frame novel view synthesis
│   ├── val_render.png              # ground truth vs NeRF render
│   └── psnr_curve.png              # training PSNR over 20k steps
├── test_images/                    # test images for classifier
├── accuracy.ipynb                  # checkpoint inspection notebook
├── labels.txt                      # 80 class names
├── taxonomy.csv                    # Khana taxonomy
└── README.md
```

## Large Files (Google Drive)
Too large for GitHub — download and place in `checkpoints/`:
- `khana_best.pt` — 334.5 MB — trained ConvNeXt-Base classifier
- `sam_vit_h_4b8939.pth` — 2.39 GB — SAM ViT-H weights
- `nerf_weights.pt` — 1.7 MB — trained NeRF coarse + fine MLPs

---

## Setup & Run

```bash
pip install torch torchvision timm segment-anything opencv-python pillow matplotlib imageio
```

### Problem 1 — Classification
```bash
# train
python codes/train.py

# predict on test images (place images in test_images/)
python codes/predict_khana.py
```

### Problem 2 — Detection (clean BEV images)
```bash
python codes/detect_thali.py path/to/thali.jpg --bev none
```

### Problem 3 — BEV + Detection (natural images)
```bash
# interactive: click 4 tray corners in popup window, press Enter
python codes/detect_thali.py path/to/thali.jpg --bev interactive

# auto: hough circle detection
python codes/detect_thali.py path/to/thali.jpg --bev auto
```

### Bonus — NeRF
```bash
# calibrate camera (need 30+ aruco grid photos)
python codes/step0_calibrate.py

# estimate poses (need 30+ thali photos with aruco marker visible)
python codes/step0_pose_and_dataset.py

# 2D neural field warmup
python codes/step1_neural_field_2d.py --image path/to/image.jpg

# full 3D NeRF (run on GPU, ~2 hours on T4)
python codes/step2_nerf_3d.py
```

### Check All Results
```bash
python codes/check_results.py
```

---

## References
- Liu et al. (2022) — [A ConvNet for the 2020s](https://arxiv.org/abs/2201.03545)
- Kirillov et al. (2023) — [Segment Anything](https://arxiv.org/abs/2304.02643)
- Mildenhall et al. (2020) — [NeRF: Representing Scenes as Neural Radiance Fields](https://arxiv.org/abs/2003.08934)
- Gaur et al. (2025) — [Khana: A Comprehensive Indian Cuisine Dataset](https://arxiv.org/pdf/2509.06006)
