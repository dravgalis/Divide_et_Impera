#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Actor.h"

#include "ProceduralMeshComponent.h"
#include "Materials/MaterialInterface.h"
#include "Engine/Texture2D.h"

#include "HexMapGenerator.generated.h"

class UProceduralMeshComponent;

UCLASS()
class DIVIDE_ET_IMPERA_API AHexMapGenerator : public AActor
{
    GENERATED_BODY()

public:
    AHexMapGenerator();

#if WITH_EDITOR
    virtual void OnConstruction(const FTransform& Transform) override;
#endif
    virtual void BeginPlay() override;

    // ===== MAP SETTINGS =====
    UPROPERTY(EditAnywhere, Category = "Hex|Map")
    int32 MapWidth = 128;

    UPROPERTY(EditAnywhere, Category = "Hex|Map")
    int32 MapHeight = 80;

    // ===== TILE SETTINGS =====
    UPROPERTY(EditAnywhere, Category = "Hex|Tile", meta = (ClampMin = "1.0"))
    float HexRadius = 100.0f; // pointy-top radius (center -> corner)

    // ВАЖНО: 1 = просто 6 треугольников; 6 = плотная сетка внутри гекса
    UPROPERTY(EditAnywhere, Category = "Hex|Tile", meta = (ClampMin = "1", ClampMax = "32"))
    int32 Subdivisions = 5;
    UPROPERTY(EditAnywhere, Category = "Hex|Tile", meta = (ClampMin = "0.0"))
    float ShoreXY = 6.0f;   // насколько второй слой уходит внутрь (uu)

    UPROPERTY(EditAnywhere, Category = "Hex|Tile", meta = (ClampMin = "0.0", ClampMax = "1.0"))
    float ShoreZAlpha = 0.5f; // насколько второй слой опускается к BoundaryZ (0..1)
    UPROPERTY(EditAnywhere, Category = "Hex|Tile", meta = (ClampMin = "0.1"))
    float WeldEpsXY = 15.0f; // uu, радиус "схлопа" внутри ОДНОГО гекса
    // ===== HEIGHTS PER LEVEL (0..8) =====
// Заполни в Editor как хочешь: например Coast ниже Plains и т.д.
    UPROPERTY(EditAnywhere, Category = "Hex|Height")
    TArray<float> LevelHeights = { 0.f, 50.f, 120.f, 180.f, 140.f, 200.f, 160.f, 320.f, 20.f };

    // ===== (FUTURE) HEIGHTMAP =====
    UPROPERTY(EditAnywhere, Category = "Hex|Map")
    FString HeightmapFilePath = TEXT("S:/Divide et Impera/Divide_et_Impera/Content/MapsData/earth_terrains_raw.txt");

protected:
    UPROPERTY(VisibleAnywhere, Category = "Hex|Components")
    TObjectPtr<UProceduralMeshComponent> ProcMesh;
    // Материал, который будет на ProcMesh (укажешь в Editor)
    UPROPERTY(EditAnywhere, Category = "Hex|Material")
    UMaterialInterface* MapMaterial = nullptr;

    // Текстура-карта высот (создаётся в рантайме)
    UPROPERTY(Transient)
    UTexture2D* HeightMapTex = nullptr;
private:
    void GenerateFlatMeshOneChunk();
    UTexture2D* BuildHeightTexture(const TArray<int32>& TileHeights, int32 W, int32 H);
    void ApplyMaterialAndTexture();
    // Axial (q,r) -> world XY (pointy-top)
    FVector2D AxialToWorld(int32 Q, int32 R) const;

    void GetHexCornersPointy(const FVector2D& Center, TArray<FVector2D>& OutCorners) const;

    // Vertex dedupe
    int32 GetOrAddVertex(
        const FVector& P,
        TArray<FVector>& Verts,
        TArray<FVector>& Normals,
        TArray<FVector2D>& UVs,
        TArray<FProcMeshTangent>& Tangents,
        TMap<FIntVector, int32>& VertexMap,
        float UVScale
    ) const;

};
