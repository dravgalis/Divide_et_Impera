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
            const FVector Center3D(Center2D.X, Center2D.Y, TileZ);
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

            // ---------- shared second-ring corner vertices (one per corner) ----------
            int32 Corner2Index[6];

            for (int32 C = 0; C < 6; ++C)
            {
                // базовая точка = предпоследняя на радиусе к углу (k = N-1)
                FVector P = Vertices[RayIndex[C][N - 1]];

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
                const int32 SideOpp = (Side + 3) % 6;

                if (GetNeighbor_EvenR_SideOrder(Q, R, SideOpp, NQ, NR))
                {
                    if (bHeightsOk && NQ >= 0 && NQ < MapWidth && NR >= 0 && NR < MapHeight)
                    {
                        const int32 NeighborLevel = TileHeights[NR * MapWidth + NQ];
                        NeighborZ = GetZFromLevel(NeighborLevel, LevelHeights);
                    }
                }

                // --- "высокий принимает Z низкого" (граница = минимум) ---
                const float BoundaryZ = FMath::Min(TileZ, NeighborZ);
                const FVector2D A2 = Corners[Side];
                const FVector2D B2 = Corners[(Side + 1) % 6];

                AddWallStrip(Center2D, A2, B2, TileZ, BoundaryZ);

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
                        // --- FIX: на внешнем ребре XY обязаны быть строго на линии ребра (общие для двух тайлов) ---
                        if (bOnOuterEdge)
                        {
                            // на ребре всегда i+j=N => j = N-i
                            const float tEdge = float(j) / float(N); // 0..1 вдоль A->B
                            const FVector2D EdgeXY = FMath::Lerp(A2, B2, tEdge);

                            P.X = EdgeXY.X;
                            P.Y = EdgeXY.Y;
                        }
                        P.Z = bOnOuterEdge ? BoundaryZ : TileZ;

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
                            Idx[i][j] = AddVertexSimple(P, TileId01, Ring01Local);
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
