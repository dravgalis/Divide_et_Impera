#include "HexMapGenerator.h"

#include "ProceduralMeshComponent.h"
#include "Kismet/KismetMathLibrary.h"
#include "Misc/FileHelper.h"
#include "Math/IntPoint.h"
#include "Materials/MaterialInstanceDynamic.h"
#include "Engine/Texture2D.h"
#include "PixelFormat.h"

AHexMapGenerator::AHexMapGenerator()
{
    PrimaryActorTick.bCanEverTick = false;

    ProcMesh = CreateDefaultSubobject<UProceduralMeshComponent>(TEXT("ProcMesh"));
    SetRootComponent(ProcMesh);

    ProcMesh->bUseAsyncCooking = true;
    ProcMesh->SetCollisionEnabled(ECollisionEnabled::QueryAndPhysics);
    ProcMesh->SetCollisionObjectType(ECC_WorldStatic);
}

#if WITH_EDITOR
void AHexMapGenerator::OnConstruction(const FTransform& Transform)
{
    Super::OnConstruction(Transform);
    GenerateFlatMeshOneChunk();
}
#endif

void AHexMapGenerator::BeginPlay()
{
    Super::BeginPlay();
    GenerateFlatMeshOneChunk();
}

// -------------------- Hex math (pointy-top) --------------------
FVector2D AHexMapGenerator::AxialToWorld(int32 Q, int32 R) const
{
    // EVEN-R offset grid (pointy-top)
    const float Width = FMath::Sqrt(3.f) * HexRadius;
    const float Height = 1.5f * HexRadius;

    const float X = Width * (Q + (R % 2 == 0 ? 0.f : 0.5f));
    const float Y = Height * R;

    // --- зеркалим ОТНОСИТЕЛЬНО ЦЕНТРА карты ---
    const float MapWorldWidth = Width * (MapWidth + 0.5f);
    const float CenterX = MapWorldWidth * 0.5f;

    const float MirroredX = CenterX * 2.f - X;

    return FVector2D(MirroredX, Y);
}

void AHexMapGenerator::GetHexCornersPointy(const FVector2D& Center, TArray<FVector2D>& OutCorners) const
{
    OutCorners.SetNum(6);

    // pointy-top: angle offset -30 degrees
    for (int32 i = 0; i < 6; ++i)
    {
        const float AngleDeg = 60.f * i - 30.f;
        const float AngleRad = FMath::DegreesToRadians(AngleDeg);
        OutCorners[i] = Center + FVector2D(FMath::Cos(AngleRad), FMath::Sin(AngleRad)) * HexRadius;
    }
}
static bool LoadHeightmapCSV(const FString& Path, int32 MapWidth, int32 MapHeight, TArray<int32>& OutHeights)
{
    OutHeights.Reset();

    FString FileData;
    if (!FFileHelper::LoadFileToString(FileData, *Path))
    {
        UE_LOG(LogTemp, Error, TEXT("Failed to load heightmap: %s"), *Path);
        return false;
    }

    FileData.ReplaceInline(TEXT("\n"), TEXT(""));
    FileData.ReplaceInline(TEXT("\r"), TEXT(""));

    TArray<FString> Tokens;
    FileData.ParseIntoArray(Tokens, TEXT(","), true);

    const int32 Expected = MapWidth * MapHeight;
    if (Tokens.Num() != Expected)
    {
        UE_LOG(LogTemp, Error, TEXT("Heightmap size mismatch. Got %d, expected %d"), Tokens.Num(), Expected);
        return false;
    }

    OutHeights.Reserve(Expected);
    for (const FString& T : Tokens)
    {
        OutHeights.Add(FCString::Atoi(*T));
    }

    UE_LOG(LogTemp, Log, TEXT("Heightmap loaded OK (%d values)"), OutHeights.Num());
    return true;
}

// -------------------- Vertex dedupe --------------------

int32 AHexMapGenerator::GetOrAddVertex(
    const FVector& P,
    TArray<FVector>& Verts,
    TArray<FVector>& Normals,
    TArray<FVector2D>& UVs,
    TArray<FProcMeshTangent>& Tangents,
    TMap<FIntVector, int32>& VertexMap,
    float UVScale
) const
{
    const float Quant = 10.f; // 0.1 uu
    const FIntVector Key(
        FMath::RoundToInt(P.X * Quant),
        FMath::RoundToInt(P.Y * Quant),
        FMath::RoundToInt(P.Z * Quant)
    );

    if (const int32* Found = VertexMap.Find(Key))
    {
        return *Found;
    }

    const int32 NewIndex = Verts.Add(P);

    Normals.Add(FVector(0, 0, 1));
    UVs.Add(FVector2D(P.X * UVScale, P.Y * UVScale));
    Tangents.Add(FProcMeshTangent(1, 0, 0));

    VertexMap.Add(Key, NewIndex);
    return NewIndex;
}

static bool GetNeighbor_EvenR_SideOrder(int32 Col, int32 Row, int32 Side, int32& OutCol, int32& OutRow)
{
    // Side: 0..5 соответствует твоему edge Corners[Side] -> Corners[Side+1]
    // 0=E, 1=NE, 2=NW, 3=W, 4=SW, 5=SE  (pointy-top, even-r)

    const bool bOdd = (Row & 1) != 0;

    // Для even-r (odd rows shifted right):
    // even row: NE=(0,-1) NW=(-1,-1) SE=(0,+1) SW=(-1,+1)
    // odd  row: NE=(+1,-1) NW=(0,-1) SE=(+1,+1) SW=(0,+1)

    switch (Side)
    {
    case 0: // E
        OutCol = Col + 1; OutRow = Row; return true;

    case 3: // W
        OutCol = Col - 1; OutRow = Row; return true;

    case 1: // NE
        OutCol = Col + (bOdd ? 1 : 0);
        OutRow = Row - 1;
        return true;

    case 2: // NW
        OutCol = Col + (bOdd ? 0 : -1);
        OutRow = Row - 1;
        return true;

    case 5: // SE
        OutCol = Col + (bOdd ? 1 : 0);
        OutRow = Row + 1;
        return true;

    case 4: // SW
        OutCol = Col + (bOdd ? 0 : -1);
        OutRow = Row + 1;
        return true;
    }

    OutCol = Col; OutRow = Row;
    return false;
}

static float GetZFromLevel(int32 Level, const TArray<float>& LevelHeights)
{
    const int32 L = FMath::Clamp(Level, 0, 8);
    if (LevelHeights.Num() >= 9)
    {
        return LevelHeights[L];
    }
    // запасной вариант, если массив не заполнен
    return float(L) * 100.f;
}
static float Hash01_U32(uint32 x)
{
    x ^= x >> 16;
    x *= 0x7feb352d;
    x ^= x >> 15;
    x *= 0x846ca68b;
    x ^= x >> 16;
    return float(x & 0x00FFFFFF) / float(0x01000000);
}

static float HashSigned_U32(uint32 x) // [-1..+1]
{
    return Hash01_U32(x) * 2.f - 1.f;
}

// общий шум по миру: одинаковый для одинаковых XY (важно для стыков между гексами)
static float SmoothStep(float t)
{
    // плавная кривая 0..1 (без резких границ)
    return t * t * (3.f - 2.f * t);
}

static float WorldNoiseSigned_Smooth(const FVector& P, float CellSize, uint32 Seed)
{
    const float Inv = (CellSize > 0.001f) ? (1.f / CellSize) : 1.f;

    const float fx = P.X * Inv;
    const float fy = P.Y * Inv;

    const int32 x0 = FMath::FloorToInt(fx);
    const int32 y0 = FMath::FloorToInt(fy);
    const int32 x1 = x0 + 1;
    const int32 y1 = y0 + 1;

    float tx = fx - float(x0);
    float ty = fy - float(y0);

    tx = SmoothStep(FMath::Clamp(tx, 0.f, 1.f));
    ty = SmoothStep(FMath::Clamp(ty, 0.f, 1.f));

    auto Sample = [&](int32 gx, int32 gy)->float
        {
            uint32 h =
                2166136261u ^
                (uint32(gx) * 16777619u) ^
                (uint32(gy) * 2166136261u) ^
                (Seed * 374761393u);

            return HashSigned_U32(h); // [-1..1]
        };

    const float v00 = Sample(x0, y0);
    const float v10 = Sample(x1, y0);
    const float v01 = Sample(x0, y1);
    const float v11 = Sample(x1, y1);

    const float a = FMath::Lerp(v00, v10, tx);
    const float b = FMath::Lerp(v01, v11, tx);
    return FMath::Lerp(a, b, ty); // [-1..1], но уже плавно
}
static float Smoothstep01(float x)
{
    x = FMath::Clamp(x, 0.f, 1.f);
    return x * x * (3.f - 2.f * x);
}

// Плавный "купол + плато":
// - plateau01 = доля радиуса (0..1), которая считается "плато" (например 0.25..0.35)
// - дальше спад по smoothstep
static float DomeWithPlateau01(float r01, float plateau01, float fallPower)
{
    r01 = FMath::Clamp(r01, 0.f, 1.f);

    if (r01 <= plateau01)
        return 1.f; // плато

    // нормализуем от plateau..1 => 0..1
    float t = (r01 - plateau01) / FMath::Max(1e-6f, (1.f - plateau01));
    t = Smoothstep01(t);

    // fallPower: >1 = круче к краю, <1 = мягче
    t = FMath::Pow(t, fallPower);
    return 1.f - t; // 1 на плато -> 0 на краю
}

static float RidgedFBM(const FVector& P, float BaseCellSize, uint32 Seed, int32 Octaves, float RidgeSharpness)
{
    Octaves = FMath::Clamp(Octaves, 1, 5);

    float sum = 0.f;
    float amp = 1.f;
    float cell = BaseCellSize;

    for (int32 o = 0; o < Octaves; ++o)
    {
        const float n = WorldNoiseSigned_Smooth(P, cell, Seed ^ (uint32)(0x9E3779B9u + o * 1013u)); // [-1..1]

        // ridged: 1 - abs(noise) => пики/гребни
        float r = 1.f - FMath::Abs(n);          // [0..1]
        r = FMath::Pow(FMath::Clamp(r, 0.f, 1.f), RidgeSharpness); // усилить "гребни"

        sum += r * amp;

        amp *= 0.5f;
        cell *= 0.5f; // больше деталей на следующих октавах
    }

    // нормируем примерно в 0..1
    const float norm = 1.f - FMath::Pow(0.5f, (float)Octaves);
    return (norm > 0.0001f) ? (sum / norm) : sum;
}

UTexture2D* AHexMapGenerator::BuildHeightTexture(const TArray<int32>& TileHeights, int32 W, int32 H)
{
    if (TileHeights.Num() != W * H) return nullptr;

    // 1 канал (байт) на пиксель
    UTexture2D* Tex = UTexture2D::CreateTransient(W, H, PF_G8);
    Tex->MipGenSettings = TMGS_NoMipmaps;
    Tex->CompressionSettings = TC_Grayscale;
    Tex->SRGB = false;
    Tex->Filter = TF_Nearest;
    Tex->AddressX = TA_Clamp;
    Tex->AddressY = TA_Clamp;

    // Заполняем mip0
    FTexture2DMipMap& Mip = Tex->GetPlatformData()->Mips[0];
    void* Data = Mip.BulkData.Lock(LOCK_READ_WRITE);

    uint8* Bytes = (uint8*)Data;
    for (int32 y = 0; y < H; ++y)
    {
        for (int32 x = 0; x < W; ++x)
        {
            int32 Level = TileHeights[y * W + x]; // 0..8
            Level = FMath::Clamp(Level, 0, 8);

            // Кодируем 0..8 в 0..255 (шаг 32)
            Bytes[y * W + x] = (uint8)(Level * 32);
        }
    }

    Mip.BulkData.Unlock();
    Tex->UpdateResource();

    return Tex;
}

void AHexMapGenerator::ApplyMaterialAndTexture()
{
    if (!ProcMesh || !MapMaterial || !HeightMapTex) return;

    UMaterialInstanceDynamic* MID = ProcMesh->CreateDynamicMaterialInstance(0, MapMaterial);
    if (!MID) return;

    // В материале должен быть Texture Parameter с именем HeightMapTex
    MID->SetTextureParameterValue(TEXT("HeightMapTex"), HeightMapTex);

    // Эти параметры тоже полезны (если будешь считать гексы в материале)
    MID->SetScalarParameterValue(TEXT("MapWidth"), (float)MapWidth);
    MID->SetScalarParameterValue(TEXT("MapHeight"), (float)MapHeight);
    MID->SetScalarParameterValue(TEXT("HexRadius"), HexRadius);

    // Pivot для твоего зеркала X (как в AxialToWorld)
    const float Width = FMath::Sqrt(3.f) * HexRadius;
    const float MapWorldWidth = Width * (MapWidth + 0.5f);
    const float CenterX = MapWorldWidth * 0.5f;
    MID->SetScalarParameterValue(TEXT("MirrorPivotX"), CenterX);
}

// -------------------- Mesh generation --------------------

void AHexMapGenerator::GenerateFlatMeshOneChunk()
{
    if (!ProcMesh) return;
    ProcMesh->ClearAllMeshSections();

    const int32 N = FMath::Max(1, Subdivisions);
    TArray<int32> TileHeights;
    const bool bHeightsOk = LoadHeightmapCSV(HeightmapFilePath, MapWidth, MapHeight, TileHeights);

    TArray<FVector> Vertices;
    TArray<int32> Triangles;
    TArray<FVector> Normals;
    TArray<FVector2D> UV0;
    TArray<FProcMeshTangent> Tangents;
    TArray<FLinearColor> VertexColors; // optional
    VertexColors.Reset();

    // Подстрой UV-скейла так, чтобы UV были адекватные по величине
    // (можешь заменить на свою логику)
    const float ApproxWorldWidth = HexRadius * FMath::Sqrt(3.f) * (MapWidth + MapHeight * 0.5f);
    const float UVScale = (ApproxWorldWidth > 1.f) ? (1.f / ApproxWorldWidth) : 0.001f;
    auto AddVertexSimple = [&](const FVector& P, float TileId01, float Border01)->int32
        {
            const int32 Idx = Vertices.Add(P);
            Normals.Add(FVector(0, 0, 1));
            UV0.Add(FVector2D(P.X * UVScale, P.Y * UVScale));
            Tangents.Add(FProcMeshTangent(1, 0, 0));

            Border01 = FMath::Clamp(Border01, 0.f, 1.f);
            VertexColors.Add(FLinearColor(TileId01, Border01, 0, 1)); // R=TileID, G=BorderBlend

            return Idx;
        };
    
    // Карта для дедупа вершин по позиции => общие границы между гексами
    TMap<FIntVector, int32> VertexMap;
    VertexMap.Reserve(MapWidth * MapHeight * 16);
    auto AddVertexSharedXY = [&](const FVector& P, float TileId01, float Border01)->int32
        {
            const float Quant = 10.f; // 0.1uu
            const FIntVector Key(
                FMath::RoundToInt(P.X * Quant),
                FMath::RoundToInt(P.Y * Quant),
                0 // Z игнорируем специально
            );

            if (int32* Found = VertexMap.Find(Key))
            {
                const int32 Idx = *Found;

                // ВАЖНО: раз вершина общая — Z должен стать тем, что мы сейчас посчитали
                Vertices[Idx].Z = P.Z;

                // Цвет: Border берём max, TileId можно усреднить
                Border01 = FMath::Clamp(Border01, 0.f, 1.f);
                if (VertexColors.IsValidIndex(Idx))
                {
                    VertexColors[Idx].G = FMath::Max(VertexColors[Idx].G, Border01);
                    if (FMath::Abs(VertexColors[Idx].R - TileId01) > 0.001f)
                        VertexColors[Idx].R = 0.5f * (VertexColors[Idx].R + TileId01);
                }

                return Idx;
            }

            const int32 NewIdx = Vertices.Add(P);

            Normals.Add(FVector(0, 0, 1));
            UV0.Add(FVector2D(P.X * UVScale, P.Y * UVScale));
            Tangents.Add(FProcMeshTangent(1, 0, 0));

            Border01 = FMath::Clamp(Border01, 0.f, 1.f);
            VertexColors.Add(FLinearColor(TileId01, Border01, 0, 1));

            VertexMap.Add(Key, NewIdx);
            return NewIdx;
        };

    // Резервы (грубая оценка)
    Vertices.Reserve(MapWidth * MapHeight * (6 * (N + 1) * (N + 2) / 2));
    Triangles.Reserve(MapWidth * MapHeight * (6 * (N * N) * 3 * 2)); // приблизительно

    TArray<FVector2D> Corners;
    Corners.Reserve(6);
    auto AddWallStrip = [&](const FVector2D& TileCenter2D, const FVector2D& A2, const FVector2D& B2, float TopZ, float BottomZ)
        {
            if (TopZ <= BottomZ) return; // стенку строит только более высокий тайл

            // Определяем "наружу" (от центра тайла)
            const FVector2D Mid2 = (A2 + B2) * 0.5f;
            const FVector2D EdgeDir = (B2 - A2).GetSafeNormal();

            // Две возможные нормали в 2D, выберем ту, что смотрит ОТ центра
            FVector2D N2(EdgeDir.Y, -EdgeDir.X);
            if (FVector2D::DotProduct(N2, Mid2 - TileCenter2D) < 0.f)
            {
                N2 = -N2;
            }
            const FVector WallNormal3(N2.X, N2.Y, 0.f);

            for (int32 k = 0; k < N; ++k)
            {
                const float t0 = float(k) / float(N);
                const float t1 = float(k + 1) / float(N);

                const FVector2D P0_2 = FMath::Lerp(A2, B2, t0);
                const FVector2D P1_2 = FMath::Lerp(A2, B2, t1);

                const FVector T0(P0_2.X, P0_2.Y, TopZ);
                const FVector T1(P1_2.X, P1_2.Y, TopZ);
                const FVector B1(P1_2.X, P1_2.Y, BottomZ);
                const FVector B0(P0_2.X, P0_2.Y, BottomZ);

                const int32 i0 = Vertices.Add(T0);
                const int32 i1 = Vertices.Add(T1);
                const int32 i2 = Vertices.Add(B1);
                const int32 i3 = Vertices.Add(B0);

                // Нормали стенки наружу
                Normals.Add(WallNormal3);
                Normals.Add(WallNormal3);
                Normals.Add(WallNormal3);
                Normals.Add(WallNormal3);

                UV0.Add(FVector2D(0, 0));
                UV0.Add(FVector2D(1, 0));
                UV0.Add(FVector2D(1, 1));
                UV0.Add(FVector2D(0, 1));

                Tangents.Add(FProcMeshTangent(1, 0, 0));
                Tangents.Add(FProcMeshTangent(1, 0, 0));
                Tangents.Add(FProcMeshTangent(1, 0, 0));
                Tangents.Add(FProcMeshTangent(1, 0, 0));

                // Треугольники так, чтобы лицевая сторона была наружу
                Triangles.Add(i0); Triangles.Add(i1); Triangles.Add(i2);
                Triangles.Add(i0); Triangles.Add(i2); Triangles.Add(i3);
            }
        };

    // Генерим сетку гексов параллелограммом по axial (q,r)
    for (int32 R = 0; R < MapHeight; ++R)
    {
        for (int32 Q = 0; Q < MapWidth; ++Q)
        {
            const FVector2D Center2D = AxialToWorld(Q, R);
            GetHexCornersPointy(Center2D, Corners);
            int32 Level = 0;
            if (bHeightsOk)
            {
                Level = TileHeights[R * MapWidth + Q]; // 0..8
            }
            const float TileId01 = float(Level) / 8.0f;
            const float TileZ = GetZFromLevel(Level, LevelHeights);
            const bool bShoreLevel = (Level == 2 || Level == 4 || Level == 6);
            // ---- стабильный seed на гекс ----
            const uint32 TileSeed =
                73856093u * uint32(Q) ^
                19349663u * uint32(R) ^
                83492791u * uint32(Level);

            const FVector Center3D(Center2D.X, Center2D.Y, TileZ);
            // ---------- Shared vertices inside ONE hex (no duplicates between sectors) ----------
            auto ComputeMountainZ = [&](const FVector& P,
                float r01,
                bool bOnOuterEdge,
                bool bNeighborIsMountain,
                float TileZLocal,
                float NeighborZLocal,
                float BoundaryZLocal,
                uint32 TileSeedLocal,
                uint32 EdgeSeedLocal) -> float
                {
                    // ---- сильно мягче профиль ----
                    const float Plateau01 = 0.42f;     // было ~0.35 -> больше плато = меньше остроты
                    const float FallPower = 0.90f;     // <1 = мягче и "шире" гора

                    const float PeakExtra = MountainPeakExtraZ;

                    // ---- шум мягче ----
                    const float BaseCell = MountainNoiseCellSize * 1.25f;  // крупнее детали
                    const int32 Oct = FMath::Clamp(MountainOctaves, 1, 4); // меньше октав = меньше "зубов"
                    const float Ridge = FMath::Max(1.10f, MountainRidgeSharpness * 0.70f); // меньше резкость

                    const float RandA = Hash01_U32(TileSeedLocal ^ 0xB5297A4Du);
                    const float RandB = Hash01_U32(TileSeedLocal ^ 0x68E31DA4u);

                    const float PeakMul = FMath::Lerp(0.95f, 1.18f, RandA);
                    const float NoiseMul = FMath::Lerp(0.55f, 0.95f, RandB);

                    // профиль (плато + плавный спад)
                    const float Profile01 = DomeWithPlateau01(r01, Plateau01, FallPower); // 1 в центре, 0 на краю

                    // ---- цель по краю ----
                    float EdgeTargetZ = BoundaryZLocal;

                    if (bNeighborIsMountain)
                    {
                        // базовая "склейка" 7↔7
                        EdgeTargetZ = 0.5f * (TileZLocal + NeighborZLocal);
                        const float Zmin = FMath::Min(TileZLocal, NeighborZLocal);
                        const float Zmax = FMath::Max(TileZLocal, NeighborZLocal);
                        EdgeTargetZ = FMath::Clamp(EdgeTargetZ, Zmin, Zmax);
                    }

                    const float PeakZ = TileZLocal + PeakExtra * PeakMul;

                    // BaseZ: от EdgeTargetZ (у края) к PeakZ (в центре)
                    float BaseZ = FMath::Lerp(EdgeTargetZ, PeakZ, Profile01);

                    // ---- шум: сильно гасим к краю и смягчаем форму ----
                    float Fade = Profile01;
                    Fade = Fade * Fade * Fade; // ещё сильнее гасим к краю

                    float NoiseAmp = MountainNoiseAmp * NoiseMul * Fade;

                    // На самом ребре шум почти выключаем (стыки)
                    if (bOnOuterEdge)
                    {
                        NoiseAmp *= (bNeighborIsMountain ? 0.10f : 0.0f);
                    }

                    const uint32 NoiseSeed = (bNeighborIsMountain ? (EdgeSeedLocal ^ 0xA24BAEDDu) : (TileSeedLocal ^ 0xA24BAEDDu));

                    const float ridged01 = RidgedFBM(P, BaseCell, NoiseSeed, Oct, Ridge); // 0..1

                    // смягчаем ridged: вместо "острых гребней" делаем мягче
                    const float rSoft = Smoothstep01(ridged01);          // 0..1
                    const float nSigned = (rSoft - 0.5f) * 2.f;          // [-1..1]

                    float Z = BaseZ + nSigned * NoiseAmp;

                    // =========================================================
                    // 7↔7 СКЛЕЙКА В ПОЛОСЕ У РЕБРА (не только на самом ребре)
                    // =========================================================
                    if (bNeighborIsMountain)
                    {
                        // ширина полосы у края (в долях радиуса)
                        const float EdgeBand01 = 0.18f; // 0.12..0.25: больше = сильнее сглаживание края

                        // edgeT = 0 внутри, 1 у края
                        float edgeT = (r01 - (1.f - EdgeBand01)) / FMath::Max(1e-6f, EdgeBand01);
                        edgeT = Smoothstep01(edgeT);

                        // Общий "идеальный" Z у края (симметричный для пары)
                        const float EdgeBase = 0.5f * (TileZLocal + NeighborZLocal);

                        // очень слабая вариация на ребре (чтобы не резать)
                        const float eRidged = RidgedFBM(P, BaseCell, EdgeSeedLocal ^ 0x77E35A1Du, FMath::Max(1, Oct - 1), 1.15f);
                        const float eSoft = Smoothstep01(eRidged);
                        const float eSigned = (eSoft - 0.5f) * 2.f;

                        const float EdgeNoiseAmp = MountainNoiseAmp * 0.08f; // маленькое!
                        const float EdgeZ = EdgeBase + eSigned * EdgeNoiseAmp;

                        // Плавно тянем вершины к EdgeZ в полосе у края
                        Z = FMath::Lerp(Z, EdgeZ, edgeT);
                    }

                    // если сосед НЕ гора — край строго boundary
                    if (!bNeighborIsMountain)
                    {
                        if (bOnOuterEdge) Z = BoundaryZLocal;
                        Z = FMath::Max(BoundaryZLocal, Z);
                    }
                    else
                    {
                        // если сосед гора — держим в разумном диапазоне
                        const float Zmin = FMath::Min(TileZLocal, NeighborZLocal);
                        const float Zmax = FMath::Max(TileZLocal, NeighborZLocal) + PeakExtra;
                        Z = FMath::Clamp(Z, Zmin, Zmax);
                    }

                    return Z;
                };

// 1) Считаем NeighborZ и BoundaryZ для всех 6 сторон (нужно для CornerZ)
            float NeighborZSide[6];
            float BoundaryZSide[6];

            for (int32 S = 0; S < 6; ++S)
            {
                int32 NQ = Q, NR = R;
                float NeighborZ = TileZ;

                const int32 SideOpp = (S + 3) % 6; // твоя правка из-за зеркала

                if (GetNeighbor_EvenR_SideOrder(Q, R, SideOpp, NQ, NR))
                {
                    if (bHeightsOk && NQ >= 0 && NQ < MapWidth && NR >= 0 && NR < MapHeight)
                    {
                        const int32 NLvl = TileHeights[NR * MapWidth + NQ];
                        NeighborZ = GetZFromLevel(NLvl, LevelHeights);
                    }
                }

                NeighborZSide[S] = NeighborZ;
                BoundaryZSide[S] = FMath::Min(TileZ, NeighborZ);
            }

            // 2) Высота углов (каждый угол принадлежит 2 сторонам)
            float CornerZ[6];
            for (int32 C = 0; C < 6; ++C)
            {
                const int32 S0 = (C + 5) % 6; // сторона "слева" от угла
                const int32 S1 = C;           // сторона "справа" от угла
                CornerZ[C] = FMath::Min(TileZ, FMath::Min(BoundaryZSide[S0], BoundaryZSide[S1]));
            }

            const int32 CenterIndex = AddVertexSimple(Center3D, TileId01, 0.0f);
            // ===== MOUNTAIN: raise center vertex =====
            if (Level == 7)
            {
                // пер-тайл рандом как у тебя в Mountain-коде
                const float RandA = Hash01_U32(TileSeed ^ 0xB5297A4Du);
                const float PeakMul = FMath::Lerp(0.80f, 1.40f, RandA);

                FVector P = Vertices[CenterIndex];

                // БАЗА: центр = TileZ + Peak
                float Zc = TileZ + MountainPeakExtraZ * PeakMul;

                // Мягкий шум для центра (НЕ ridged), чтобы не было иглы
                const float n = WorldNoiseSigned_Smooth(P, MountainNoiseCellSize * 2.0f, TileSeed ^ 0x51ED270Bu);
                Zc += n * (MountainNoiseAmp * 0.25f); // 25% амплитуды в центре — обычно достаточно

                P.Z = Zc;
                Vertices[CenterIndex] = P;
            }

            // 4) Создаём 6 радиусов: RayIndex[corner][k], k=0..N
            TArray<TArray<int32>> RayIndex;
            RayIndex.SetNum(6);

            for (int32 C = 0; C < 6; ++C)
            {
                RayIndex[C].SetNum(N + 1);
                RayIndex[C][0] = CenterIndex;

                const FVector Corner3D(Corners[C].X, Corners[C].Y, CornerZ[C]);

                for (int32 k = 1; k <= N; ++k)
                {
                    const float t = float(k) / float(N);

                    FVector P = FMath::Lerp(Center3D, Corner3D, t);

                    // ВАЖНО: чтобы не делать "скат" по радиусу — держим Z=TileZ на k=1..N-1
                    if (k < N) P.Z = TileZ;

                    const float Border01 = float(k) / float(N); // 0 в центре, 1 на границе
                    RayIndex[C][k] = AddVertexSimple(P, TileId01, Border01);
                }
            }
            if (Level == 7 && N >= 2)
            {
                for (int32 C = 0; C < 6; ++C)
                {
                    const int32 VIdx = RayIndex[C][1];
                    FVector P = Vertices[VIdx];

                    // чуть ниже центра, чтобы не было иглы
                    const float Ring01 = 1.f / float(N);
                    const float Fall = FMath::Pow(Ring01, MountainFalloffPower);

                    const float RandA = Hash01_U32(TileSeed ^ 0xB5297A4Du);
                    const float PeakMul = FMath::Lerp(0.80f, 1.40f, RandA);

                    const float CenterTargetZ = TileZ + MountainPeakExtraZ * PeakMul;
                    const float BaseZ = FMath::Lerp(CenterTargetZ, TileZ, Fall);

                    const float n = WorldNoiseSigned_Smooth(P, MountainNoiseCellSize * 2.0f, TileSeed ^ 0x51ED270Bu);
                    P.Z = BaseZ + n * (MountainNoiseAmp * 0.25f);

                    Vertices[VIdx] = P;
                }
            }

            // ===== APPLY HILL/MOUNTAIN Z TO RAY VERTICES (sector borders inside hex) =====
            const bool bRayIsHill = (Level == 3 || Level == 5);
            const bool bRayIsMountain = (Level == 7);

            if (bRayIsMountain)
            {
                // Центр тоже должен быть горой
                {
                    FVector Pc = Vertices[CenterIndex];
                    Pc.Z = TileZ + MountainPeakExtraZ; // базово поднимем, дальше можно шумом
                    Vertices[CenterIndex] = Pc;
                }

                // Лучи: k=1..N-1 (внутренние границы секторов)
                for (int32 C = 0; C < 6; ++C)
                {
                    for (int32 k = 1; k < N; ++k)
                    {
                        const int32 VIdx = RayIndex[C][k];
                        if (!Vertices.IsValidIndex(VIdx)) continue;

                        FVector P = Vertices[VIdx];

                        const float r01 = float(k) / float(N); // 0..1

                        // Это НЕ внешний край гекса, это внутренний "луч",
                        // поэтому считаем как "не outer edge".
                        const bool bOnOuterEdge = false;

                        // сосед для лучей не нужен
                        const bool bNeighborIsMountain = false;

                        const float Zm = ComputeMountainZ(
                            P,
                            r01,
                            bOnOuterEdge,
                            bNeighborIsMountain,
                            TileZ, TileZ, TileZ, // NeighborZ/BoundaryZ не важны здесь
                            TileSeed,
                            0u
                        );

                        P.Z = Zm;

                        Vertices[VIdx] = P;
                    }
                }
            }


            // ---------- shared second-ring corner vertices (one per corner) ----------
            int32 Corner2Index[6];

            for (int32 C = 0; C < 6; ++C)
            {
                // базовая точка = предпоследняя на радиусе к углу (k = N-1)
                FVector P = Vertices[RayIndex[C][N - 1]];

                // если ХОТЯ БЫ ОДНА из двух сторон у этого угла имеет соседа ниже — топим
                const int32 S0 = (C + 5) % 6;
                const int32 S1 = C;

                const bool bLower0 = (NeighborZSide[S0] < TileZ - 0.1f);
                const bool bLower1 = (NeighborZSide[S1] < TileZ - 0.1f);
                const bool bShouldInset = bShoreLevel && (bLower0 || bLower1);

                if (bShouldInset)
                {
                    FVector2D ToCenter = (Center2D - FVector2D(P.X, P.Y)).GetSafeNormal();
                    const float CornerInset = ShoreXY * 1.55f; // подбирай 0.25..0.65

                    P.X += ToCenter.X * CornerInset;
                    P.Y += ToCenter.Y * CornerInset;

                    // Z: тянем вниз к более низкой из двух сторон угла (не ниже Boundary)
                    const float CornerBoundary = FMath::Min(BoundaryZSide[S0], BoundaryZSide[S1]);
                    P.Z = FMath::Max(CornerBoundary, FMath::Lerp(TileZ, CornerBoundary, ShoreZAlpha));
                }

                const float Border01 = float(N - 1) / float(N); // второе кольцо почти у края
                Corner2Index[C] = AddVertexSimple(P, TileId01, Border01);
            }
            // second ring side vertices: [Side][k] where k = 1..N-2
            TArray<TArray<int32>> SecondRingSide;
            SecondRingSide.SetNum(6);
            for (int32 S = 0; S < 6; ++S)
            {
                SecondRingSide[S].SetNum(N - 1);
                for (int32 k = 0; k < N - 1; ++k)
                    SecondRingSide[S][k] = INDEX_NONE;
            }

            // Каждый гекс: 6 треугольных секторов (Center, Corner[i], Corner[i+1])
            for (int32 Side = 0; Side < 6; ++Side)
            {
                // --- найти соседа по этой стороне и вычислить его Z ---
                int32 NQ = Q, NR = R;
                float NeighborZ = TileZ; // если соседа нет — используем свой Z
                int32 NeighborLevel = Level; // по умолчанию считаем, что "такой же"

                const int32 SideOpp = (Side + 3) % 6;

                if (GetNeighbor_EvenR_SideOrder(Q, R, SideOpp, NQ, NR))
                {
                    if (bHeightsOk && NQ >= 0 && NQ < MapWidth && NR >= 0 && NR < MapHeight)
                    {
                        NeighborLevel = TileHeights[NR * MapWidth + NQ];
                        NeighborZ = GetZFromLevel(NeighborLevel, LevelHeights);
                    }
                }

                const bool bThisIsHill = (Level == 3 || Level == 5);
                const bool bNeighborIsHill = (NeighborLevel == 3 || NeighborLevel == 5);
                const bool bThisIsMountain = (Level == 7);
                const bool bNeighborIsMountain = (NeighborLevel == 7);

                // --- Symmetric Edge Seed (одинаковый для пары тайлов) ---
                const int32 ACol = Q;
                const int32 ARow = R;
                const int32 BCol = NQ;
                const int32 BRow = NR;

                // сортируем пару координат, чтобы порядок не влиял
                const bool bSwap = (BRow < ARow) || (BRow == ARow && BCol < ACol);

                const int32 MinC = bSwap ? BCol : ACol;
                const int32 MinR = bSwap ? BRow : ARow;
                const int32 MaxC = bSwap ? ACol : BCol;
                const int32 MaxR = bSwap ? ARow : BRow;

                uint32 EdgeSeed =
                    1469598103u ^
                    (uint32(MinC) * 16777619u) ^
                    (uint32(MinR) * 2166136261u) ^
                    (uint32(MaxC) * 374761393u) ^
                    (uint32(MaxR) * 668265263u) ^
                    (uint32(Side) * 2246822519u);

                // --- "высокий принимает Z низкого" (граница = минимум) ---
                const float BoundaryZ = FMath::Min(TileZ, NeighborZ);
                const bool bLowerNeighbor = (NeighborZ < TileZ - 0.1f);

                const FVector2D A2 = Corners[Side];
                const FVector2D B2 = Corners[(Side + 1) % 6];
                // Нормаль стороны, направленная ВНУТРЬ гекса (для смещения 2-го слоя)
                FVector2D EdgeDir = (B2 - A2).GetSafeNormal();
                FVector2D InN(EdgeDir.Y, -EdgeDir.X); // перпендикуляр

                const FVector2D Mid2 = (A2 + B2) * 0.5f;
                if (FVector2D::DotProduct(InN, Center2D - Mid2) < 0.f)
                {
                    InN = -InN;
                }

                const FVector A3(A2.X, A2.Y, 0.f);
                const FVector B3(B2.X, B2.Y, 0.f);

                // Индексы вершин треугольной решётки (i row, j col) внутри сектора
                TArray<TArray<int32>> Idx;
                Idx.SetNum(N + 1);

                for (int32 i = 0; i <= N; ++i)
                {
                    const int32 RowLen = (N - i) + 1;
                    Idx[i].SetNum(RowLen);

                    for (int32 j = 0; j < RowLen; ++j)
                    {
                        // Барицентрическая параметризация:
                        // P = C + (A-C)*(i/N) + (B-C)*(j/N)
                        const float fi = float(i) / float(N);
                        const float fj = float(j) / float(N);

                        FVector P =
                            Center3D +
                            (A3 - Center3D) * fi +
                            (B3 - Center3D) * fj;
                        
                        const bool bOnOuterEdge = (i + j == N);       // ребро A-B
                        const bool bSecondRing = (i + j == N - 1);   // второй слой рядом с ребром
                        const bool bIsCornerA = bOnOuterEdge && (j == 0); // A
                        const bool bIsCornerB = bOnOuterEdge && (i == 0); // B
                        const FVector2D BasePos(P.X, P.Y);
                        float Z = bOnOuterEdge ? BoundaryZ : TileZ;
                        // --- FIX: на внешнем ребре XY обязаны быть строго на линии ребра (общие для двух тайлов) ---
                        if (bOnOuterEdge)
                        {
                            // на ребре всегда i+j=N => j = N-i
                            const float tEdge = float(j) / float(N); // 0..1 вдоль A->B
                            const FVector2D EdgeXY = FMath::Lerp(A2, B2, tEdge);

                            P.X = EdgeXY.X;
                            P.Y = EdgeXY.Y;
                        }

                        // --- углы границы учитывают два ребра ---
                        if (bIsCornerA || bIsCornerB)
                        {
                            float Zmin = FMath::Min(TileZ, BoundaryZ);

                            const int32 Side2 = bIsCornerA ? ((Side + 5) % 6) : ((Side + 1) % 6);
                            const int32 Side2Opp = (Side2 + 3) % 6;

                            int32 NQ2 = Q, NR2 = R;
                            float NeighborZ2 = TileZ;

                            if (GetNeighbor_EvenR_SideOrder(Q, R, Side2Opp, NQ2, NR2))
                            {
                                if (bHeightsOk && NQ2 >= 0 && NQ2 < MapWidth && NR2 >= 0 && NR2 < MapHeight)
                                {
                                    const int32 L2 = TileHeights[NR2 * MapWidth + NQ2];
                                    NeighborZ2 = GetZFromLevel(L2, LevelHeights);
                                }
                            }

                            const float BoundaryZ2 = FMath::Min(TileZ, NeighborZ2);
                            Zmin = FMath::Min(Zmin, BoundaryZ2);

                            Z = Zmin;
                        }

                        // --- берег: только тайл 2, только второй слой, только если сосед ниже ---
                        if (bShoreLevel && bSecondRing && bLowerNeighbor)
                        {
                            // ---- утопление углов второго слоя ----
                            if (i == 0 || j == 0)
                            {
                                FVector2D ToCenter = Center2D - FVector2D(P.X, P.Y);
                                ToCenter.Normalize();

                                const float CornerInset = ShoreXY * 0.50f; // 0.35–0.55 подбирается

                                P.X += ToCenter.X * CornerInset;
                                P.Y += ToCenter.Y * CornerInset;
                            }

                            // ---- параметр вдоль стороны (0..1) ----
                            const float denom = float(N - i);
                            const float t = (denom > 0.f) ? (float(j) / denom) : 0.5f;
                            const bool bCornerA = (j == 0);
                            const bool bCornerB = (i == 0);
                            const bool bCorner = bCornerA || bCornerB;

                            // вес для середины стороны
                            const float EdgeWeight = 4.f * t * (1.f - t);

                            // вес для углов (НЕ НОЛЬ!)
                            const float CornerWeight = 0.35f;

                            // итоговый вес
                            const float Weight = bCorner ? CornerWeight : EdgeWeight;

                            // ---- hash noise (stable, tile-dependent) ----
                            auto Hash01 = [](uint32 x)
                                {
                                    x ^= x >> 16;
                                    x *= 0x7feb352d;
                                    x ^= x >> 15;
                                    x *= 0x846ca68b;
                                    x ^= x >> 16;
                                    return float(x & 0x00FFFFFF) / float(0x01000000);
                                };

                            uint32 NoiseSeed =
                                TileSeed ^
                                (uint32(Side) * 1013904223u) ^
                                (uint32(i * 17 + j * 31) * 1664525u);

                            float Noise = Hash01(NoiseSeed) * 2.f - 1.f;
                            
                            const float LocalShoreMul =
                                0.75f + Hash01(TileSeed ^ 0xA341316Cu) * 0.6f;

                            // ---- направление смещения ----
                            FVector2D MoveDir = InN;
                            float ShiftMul = 1.f;


                            if (bCornerA || bCornerB)
                            {
                                const int32 OtherSide = bCornerA ? (Side + 5) % 6 : (Side + 1) % 6;

                                FVector2D A2b = Corners[OtherSide];
                                FVector2D B2b = Corners[(OtherSide + 1) % 6];
                                FVector2D EdgeDir2 = (B2b - A2b).GetSafeNormal();
                                FVector2D InN2(EdgeDir2.Y, -EdgeDir2.X);

                                FVector2D Mid2b = (A2b + B2b) * 0.5f;
                                if (FVector2D::DotProduct(InN2, Center2D - Mid2b) < 0.f)
                                    InN2 = -InN2;

                                MoveDir = (InN + InN2).GetSafeNormal();
                                ShiftMul = 0.6f; // углы двигаем слабее
                            }

                            
                            const float Angle =
                                (Hash01(TileSeed ^ (Side * 92837111u)) - 0.5f) * 0.35f; // ±20°

                            const float c = FMath::Cos(Angle);
                            const float s = FMath::Sin(Angle);

                            MoveDir = FVector2D(
                                MoveDir.X * c - MoveDir.Y * s,
                                MoveDir.X * s + MoveDir.Y * c
                            );
                            const float BreakChance = Hash01(TileSeed ^ (Side * 912931u));
                            const float BreakMul = (BreakChance < 0.25f) ? 0.15f : 1.0f;
                            // ---- базовый shift (как было) ----
                            float FinalShift =
                                ShoreXY *
                                LocalShoreMul *
                                Weight *
                                BreakMul *
                                (1.f + Noise * 0.45f);

                            // ---- сглаживание ТОЛЬКО для вершин стороны второго кольца (НЕ углы) ----
                            if (bSecondRing && !bCorner && i > 0 && j > 0)
                            {
                                // Noise соседей вдоль ребра: (i-1, j+1) и (i+1, j-1)
                                const uint32 NoiseSeedPrev =
                                    TileSeed ^
                                    (uint32(Side) * 1013904223u) ^
                                    (uint32((i - 1) * 17 + (j + 1) * 31) * 1664525u);

                                const uint32 NoiseSeedNext =
                                    TileSeed ^
                                    (uint32(Side) * 1013904223u) ^
                                    (uint32((i + 1) * 17 + (j - 1) * 31) * 1664525u);

                                const float NoisePrev = Hash01(NoiseSeedPrev) * 2.f - 1.f;
                                const float NoiseNext = Hash01(NoiseSeedNext) * 2.f - 1.f;
                                // --- Weight для PREV (i-1, j+1) ---
                                const int32 iPrev = i - 1;
                                const int32 jPrev = j + 1;
                                const float denomPrev = float(N - iPrev);
                                const float tPrev = (denomPrev > 0.f) ? (float(jPrev) / denomPrev) : 0.5f;
                                const float WeightPrev = 4.f * tPrev * (1.f - tPrev);

                                // --- Weight для NEXT (i+1, j-1) ---
                                const int32 iNext = i + 1;
                                const int32 jNext = j - 1;
                                const float denomNext = float(N - iNext);
                                const float tNext = (denomNext > 0.f) ? (float(jNext) / denomNext) : 0.5f;
                                const float WeightNext = 4.f * tNext * (1.f - tNext);

                                const float ShiftPrev =
                                    ShoreXY * LocalShoreMul * WeightPrev * BreakMul * (1.f + NoisePrev * 0.45f);

                                const float ShiftNext =
                                    ShoreXY * LocalShoreMul * WeightNext * BreakMul * (1.f + NoiseNext * 0.45f);

                                const float ShiftAvg = (ShiftPrev + FinalShift + ShiftNext) / 3.f;

                                const float SmoothAlpha = 0.0f; // 0..1 (увеличивай — будет плавнее, но НЕ глубже)
                                const float ShiftSmoothed = FMath::Lerp(FinalShift, ShiftAvg, SmoothAlpha);

                                // Ограничиваем СКОРОСТЬ изменения, а не "зажимаем в диапазон"
                                // так сглаживание реально влияет, но не даёт резко вгонять внутрь
                                const float MaxDelta = ShoreXY * 10.35f; // подбирай: 0.15..0.6
                                FinalShift = FinalShift + FMath::Clamp(ShiftSmoothed - FinalShift, -MaxDelta, MaxDelta);

                            }

                            P.X += MoveDir.X * FinalShift;
                            P.Y += MoveDir.Y * FinalShift;

                            // ---- Z ----
                            const float TargetZ =
                                FMath::Lerp(
                                    TileZ,
                                    BoundaryZ,
                                    ShoreZAlpha * Weight * BreakMul
                                );

                            Z = FMath::Max(
                                BoundaryZ,
                                FMath::Lerp(TileZ, TargetZ, 1.f + Noise * 0.2f)
                            );

                        }
                        // ---------------- HILLS (Level 3 / 5) ----------------
                        if (bThisIsHill)
                        {
                            // --- Per-tile random (чтобы каждый холм отличался) ---
                            const float Rand01_A = Hash01_U32(TileSeed ^ 0xA341316Cu); // 0..1
                            const float Rand01_B = Hash01_U32(TileSeed ^ 0xC8013EA4u);
                            const float Rand01_C = Hash01_U32(TileSeed ^ 0xAD90777Du);

                            // множители параметров на конкретный холм
                            const float PeakMul = FMath::Lerp(0.70f, 1.35f, Rand01_A); // высота пика
                            const float NoiseMul = FMath::Lerp(0.60f, 1.50f, Rand01_B); // сила шума
                            const float CellMul = FMath::Lerp(0.75f, 1.40f, Rand01_C); // размер "клетки" шума

                            // немного рандомизируем спад: у одних холм круче, у других плавнее
                            const float FalloffMul = FMath::Lerp(0.85f, 1.25f, Hash01_U32(TileSeed ^ 0x7F4A7C15u));

                            // 0 в центре, 1 на границе
                            const float Ring01 = float(i + j) / float(N);

                            // спад (чем ближе к краю — тем сильнее)
                            const float FallPower = HillFalloffPower * FalloffMul;
                            const float Fall = FMath::Pow(FMath::Clamp(Ring01, 0.f, 1.f), FallPower);

                            // Цель по краю:
                            // - если сосед не холм: оставляем стандартный BoundaryZ (min)
                            // - если сосед холм и флаг включён: делаем мягкий стык между двумя холмами
                            float EdgeTargetZ = BoundaryZ;

                            if (bHillBlendEdgesWithHillNeighbors && bNeighborIsHill)
                            {
                                EdgeTargetZ = 0.5f * (TileZ + NeighborZ);
                                const float Zmin = FMath::Min(TileZ, NeighborZ);
                                const float Zmax = FMath::Max(TileZ, NeighborZ);
                                EdgeTargetZ = FMath::Clamp(EdgeTargetZ, Zmin, Zmax);
                            }

                            const float CenterTargetZ = TileZ + HillPeakExtraZ * PeakMul;
                            float BaseZ = FMath::Lerp(CenterTargetZ, EdgeTargetZ, Fall);

                            // шум затухает к краю (в центре больше)
                            const float Noise01 = 1.f - FMath::Clamp(Ring01, 0.f, 1.f);
                            const float NoiseFade = Noise01 * Noise01;

                            const float LocalCellSize = HillNoiseCellSize * CellMul;
                            const float WorldNoise = WorldNoiseSigned_Smooth(P, LocalCellSize, TileSeed ^ uint32(Level * 912931u));

                            float NoiseAmp = HillNoiseAmp * NoiseMul * NoiseFade;

                            // запрет трогать outer edge, кроме случая сосед тоже холм
                            if (bOnOuterEdge)
                            {
                                if (bHillBlendEdgesWithHillNeighbors && bNeighborIsHill)
                                {
                                    NoiseAmp = HillNoiseAmp * NoiseMul * HillEdgeNoiseAmp01;
                                }
                                else
                                {
                                    NoiseAmp = 0.f;
                                }
                            }

                            float Zhill = BaseZ + WorldNoise * NoiseAmp;

                            if (!(bHillBlendEdgesWithHillNeighbors && bNeighborIsHill))
                            {
                                // если сосед не холм — граница строго BoundaryZ
                                if (bOnOuterEdge) Zhill = BoundaryZ;

                                // не провалиться ниже boundary
                                Zhill = FMath::Max(BoundaryZ, Zhill);
                            }
                            else
                            {
                                // если сосед холм — держимся в разумном диапазоне
                                const float Zmin = FMath::Min(TileZ, NeighborZ);
                                const float Zmax = FMath::Max(TileZ, NeighborZ) + HillPeakExtraZ;
                                Zhill = FMath::Clamp(Zhill, Zmin, Zmax);
                            }

                            Z = Zhill;
                        }
                        else if (bThisIsMountain)
                        {
                            // r01 для вершины (0 центр, 1 край)
                            const float r01 = float(i + j) / float(N);

                            float Zm = ComputeMountainZ(
                                P,
                                r01,
                                bOnOuterEdge,
                                bNeighborIsMountain,
                                TileZ,
                                NeighborZ,
                                BoundaryZ,
                                TileSeed,
                                EdgeSeed
                            );

                            Z = Zm;
                        }

                        // ------------------------------------------------------
                        // (твой shore-код выше НЕ трогаем)
                        // ------------------------------------------------------
                        else if (!bOnOuterEdge)
                        {
                            Z = TileZ;
                        }

                        P.Z = Z;

                        const bool bCornerSecond = bSecondRing && (i == 0 || j == 0);

                        if (bCornerSecond)
                        {
                            // j==0 => угол Side, i==0 => угол Side+1
                            const int32 CornerId = (j == 0) ? Side : ((Side + 1) % 6);
                            Idx[i][j] = Corner2Index[CornerId];
                        }
                        else if (j == 0 && !bSecondRing)
                        {
                            Idx[i][j] = RayIndex[Side][i];
                        }
                        else if (i == 0 && !bSecondRing)
                        {
                            Idx[i][j] = RayIndex[(Side + 1) % 6][j];
                        }
                        else if (bSecondRing && (i > 0 && j > 0))
                        {
                            // индекс вдоль стороны (0..N-2)
                            const int32 k = j - 1;

                            if (SecondRingSide[Side][k] == INDEX_NONE)
                            {
                                const float Border01 = float(N - 1) / float(N);
                                SecondRingSide[Side][k] = AddVertexSimple(P, TileId01, Border01);
                            }

                            Idx[i][j] = SecondRingSide[Side][k];
                        }
                        else
                        {
                            const float Ring01Local = float(i + j) / float(N);

                            // Только для outer edge холм↔холм делаем общую вершину и общий Z по ребру
                            if (bThisIsHill && bOnOuterEdge && bHillBlendEdgesWithHillNeighbors && bNeighborIsHill)
                            {
                                // Симметричный Z на ребре (одинаковый у обоих тайлов)
                                const float EdgeNoise = WorldNoiseSigned_Smooth(P, HillNoiseCellSize, EdgeSeed); // [-1..1]
                                const float EdgeBaseZ = 0.5f * (TileZ + NeighborZ);
                                P.Z = EdgeBaseZ + EdgeNoise * (HillNoiseAmp * HillEdgeNoiseAmp01);

                                Idx[i][j] = AddVertexSharedXY(P, TileId01, Ring01Local);
                            }
                            else
                            {
                                Idx[i][j] = AddVertexSimple(P, TileId01, Ring01Local);
                            }
                            // общий Z и общий индекс на ребре для 7↔7
                            if (bThisIsMountain && bOnOuterEdge && bNeighborIsMountain)
                            {
                                // ВАЖНО: на ребре используем EdgeSeed (симметричный), чтобы оба тайла дали один и тот же Z
                                const float eRidged01 = RidgedFBM(P, MountainNoiseCellSize, EdgeSeed ^ 0xA24BAEDDu, MountainOctaves, MountainRidgeSharpness);
                                const float eSigned = (eRidged01 * 2.f - 1.f);

                                const float EdgeBaseZ = 0.5f * (TileZ + NeighborZ);
                                P.Z = EdgeBaseZ + eSigned * (MountainNoiseAmp * MountainEdgeNoiseAmp01);

                                Idx[i][j] = AddVertexSharedXY(P, TileId01, Ring01Local);
                            }
                            else
                            {
                                // оставляем как было (AddVertexSimple или hill↔hill и т.п.)
                            }

                        }

                    }

                }

                // Триангуляция треугольной решётки
                // Для каждого "квадратика" в треугольной сетке — 1 или 2 треугольника
                for (int32 i = 0; i < N; ++i)
                {
                    const int32 RowLen = (N - i) + 1;       // Idx[i].Num()
                    const int32 NextRowLen = (N - (i + 1)) + 1; // Idx[i+1].Num()

                    for (int32 j = 0; j < RowLen - 1; ++j)
                    {
                        const int32 V0 = Idx[i][j];
                        const int32 V1 = Idx[i + 1][j];
                        const int32 V2 = Idx[i][j + 1];

                        // CCW для нормалей вверх (+Z)
                        Triangles.Add(V0);
                        Triangles.Add(V2);
                        Triangles.Add(V1);

                        if (j + 1 < NextRowLen)
                        {
                            const int32 V3 = Idx[i + 1][j + 1];

                            Triangles.Add(V1);
                            Triangles.Add(V2);
                            Triangles.Add(V3);
                        }

                    }
                }
            }

        }
    }
    UE_LOG(LogTemp, Warning, TEXT("Verts=%d Colors=%d Tris=%d"),
        Vertices.Num(), VertexColors.Num(), Triangles.Num());

    ProcMesh->CreateMeshSection_LinearColor(
        0,
        Vertices,
        Triangles,
        Normals,
        UV0,
        VertexColors,
        Tangents,
        true // collision
    );
    // Создаём текстуру-карту высот и отдаём в материал
    HeightMapTex = BuildHeightTexture(TileHeights, MapWidth, MapHeight);
    ApplyMaterialAndTexture();

    ProcMesh->SetMeshSectionVisible(0, true);
}
