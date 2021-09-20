#include "SegmentationAlgorithms.h"

using namespace ogx;
using namespace ogx::Data;

/////////////////////// Custom functions ///////////////////////////////////////////////////////////////

/** \struct
\brief function

Function enables to check if two 3D vectors are pointing in the same direction

\return boolean: true - same direction, false - different directions
*/

bool is_same_direction(Math::Vector3D vector1, Math::Vector3D vector2)
{
	bool is_similarly_oriented = true;

	Real angle = Math::CalcAngleBetweenTwoVectors(vector1, vector2); // angle in degrees

	if (angle > (PI / 2))
		is_similarly_oriented = false;

	return is_similarly_oriented;
}

/** \struct
\brief function for std::sort



*/

bool sortbysecdesc(std::pair<int, int>& a, std::pair<int, int>& b)
{
	return a.second > b.second;
}


/////////////////////// Pre-segmentation methods ///////////////////////////////////////////////////////

/** \struct
\brief LPFC parameter calculation method

\Calculates the Local Plane Fitting Coefficient (LPFC) parameter for every point in the given point cloud

*/

struct Local_plane_fitting_coefficient : public ogx::Plugin::EasyMethod
{
	// parameters
	Data::ResourceID m_node_id; /*node used in segmentation*/
	Real LPFC_KNN_size; /*Number of nearest points in a point's neighbourhood*/

	// constructor
	Local_plane_fitting_coefficient() : EasyMethod(L"Micha³ Kossakowski", L"Calculates the Local Plane Fitting Coefficient (LPFC) for every point in the point cloud.")
	{
	}

	// add input/output parameters
	virtual void DefineParameters(ParameterBank& bank)
	{
		bank.Add(L"node_id", m_node_id, L"ID of node containing point cloud(s)").AsNode();
		bank.Add(L"LPFC KNN size", LPFC_KNN_size = 10, L"Number of K nearest points for LPFC parameter calculation");

	}

	virtual void Run(Context& context)
	{
		auto subtree = context.Project().TransTreeFindNode(m_node_id);
		// report error if given node was not found, this will stop execution of algorithm
		if (!subtree) ReportError(L"Node not found");


		// run with number of threads available on current machine, optional
		auto const thread_count = std::thread::hardware_concurrency();

		// perform calculations for each cloud in given subtree
		Clouds::ForEachCloud(*subtree, [&](Data::Clouds::ICloud& cloud, Data::Nodes::ITransTreeNode& node)
			{
				//access points in the cloud
				Data::Clouds::PointsRange all_points;
				cloud.GetAccess().GetAllPoints(all_points);

				std::vector<Data::Clouds::Point3D> points_xyz;
				all_points.GetXYZ(points_xyz);

				auto LPFCLayer = cloud.CreateLayer(L"LPFC", 0.0);
				std::vector<float> LPFCValues;

				float max_LPFC = 0;
				float sum = 0;

				// Calculate LPFC for every point in the point cloud
				for (auto i = 0; i < all_points.size(); i++)
				{
					Math::Point3D point = points_xyz[i].cast<Real>();

					// create local neighbourhood around point using KNNSearchKernel
					Data::Clouds::KNNSearchKernel search(point, LPFC_KNN_size);
					Data::Clouds::PointsRange neighbourhood;
					cloud.GetAccess().FindPoints(search, neighbourhood);

					std::vector<Data::Clouds::Point3D> neighbours_xyz;
					neighbourhood.GetXYZ(neighbours_xyz);

					auto plane = Math::CalcBestPlane3D(neighbours_xyz.begin(), neighbours_xyz.end());

					// Calculate the sum of every point's distance to the best fitting plane
					Real distances_sum = 0;
					for (auto neighbour : neighbours_xyz)
					{
						Math::Point3D neighbour_real = neighbour.cast<Real>();

						Real distance = Math::CalcPointToPlaneDistance3D(neighbour_real, plane, false);

						distances_sum += distance;
					}

					// calculate LPFC
					float LPFC = distances_sum / neighbours_xyz.size();

					LPFCValues.push_back(LPFC);
				}

				all_points.SetLayerVals(LPFCValues, *LPFCLayer);

			}, thread_count); // run with given number of threads, optional parameter, if not given will run in current thread
	}
};

/////////////////////// Segmentation algorithms ////////////////////////////////////////////////////////

/** \struct
\brief Region growing segmentation method

\Divides point cloud into segments of similarly oriented surfaces using region growing

*/

struct SegmentationByRegionGrowingWithLPFC : public ogx::Plugin::EasyMethod
{
	// parameters
	Data::ResourceID m_node_id;
	Real KNN_size;
	Real number_of_points;
	Real threshold_angle;
	Real threshold_growth_factor;
	Real LPFC_KNN_size;
	String LPFC_layer_name;
	Real min_LPFC;
	Real SR_size_threshold;
	Real SR_threshold_angle;

	// constructor
	SegmentationByRegionGrowingWithLPFC() : EasyMethod(L"Micha³ Kossakowski", L"Segmentation by a region growing method which creates regions of similarly oriented surfaces.")
	{
	}

	// add input/output parameters
	virtual void DefineParameters(ParameterBank& bank)
	{
		bank.Add(L"node_id", m_node_id, L"ID of node containing point cloud(s)").AsNode();
		bank.Add(L"KNN_size", KNN_size = 10, L"Number of K nearest points for region growing");
		bank.Add(L"number_of_points", number_of_points = 5000, L"Number of initial seed points");
		bank.Add(L"threshold_angle", threshold_angle = 15, L"Acceptable deviation of surface orientation [degrees]"); // angle in degrees
		bank.Add(L"threshold_growth_factor", threshold_growth_factor = 10, L"Minimal growth factor value for caluclating new region normal vector in region growing process [%]"); // growth factor in percent [%]
		bank.Add(L"LPFC_KNN_size", LPFC_KNN_size = 10, L"Number of K nearest points for LPFC parameter calculation");
		bank.Add(L"LPFC_layer_name", LPFC_layer_name = L"LPFC", L"Name of data layer containing LPFC parameter values");
		bank.Add(L"min_LPFC", min_LPFC = 0.1, L"Minimal LPFC parameter value for choosing initial seed points"); // minimal LPFC parameter value for selecting seed points (in mm)
		bank.Add(L"SR_size_threshold", SR_size_threshold = 300, L"Threshold (max) size value for smaller regions (SR) in number of points");
		bank.Add(L"SR_threshold_angle", SR_threshold_angle = 20, L"Threshold (max) angle value for grouping smaller regions (SR) [degrees]");

	}

	virtual void Run(Context& context)
	{
		auto subtree = context.Project().TransTreeFindNode(m_node_id);
		// report error if given node was not found, this will stop execution of algorithm
		if (!subtree) ReportError(L"Node not found");


		// run with number of threads available on current machine, optional
		auto const thread_count = std::thread::hardware_concurrency();

		// perform calculations for each cloud in given subtree
		Clouds::ForEachCloud(*subtree, [&](Data::Clouds::ICloud& cloud, Data::Nodes::ITransTreeNode& node)
			{
				auto start = std::chrono::steady_clock::now();

				ogx::Execution::Parameters params;
				params.Clear();
				params.AddParameter(L"node_id", ogx::Execution::ParameterType::PARAM_RESOURCE_ID, node.GetID());
				params.AddParameter(L"LPFC KNN size", ogx::Execution::ParameterType::PARAM_REAL, LPFC_KNN_size);

				if (context.Execution().ExecuteAlgorithmSync(L"Example.Local_plane_fitting_coefficient", params))
					OGX_LINE.Msg(User, L"Found Algorithm");
				else
					OGX_LINE.Msg(User, L"ERROR: Cant find the algorithm Example.Local_plane_fitting_coefficient");

				auto end = std::chrono::steady_clock::now();
				float elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
				OGX_LINE.Format(Info, L"Elapsed time LPFC: %f ms", elapsed_time);

				//access points in the cloud
				Data::Clouds::PointsRange all_points;
				cloud.GetAccess().GetAllPoints(all_points);

				auto vector_of_layers = cloud.FindLayers(LPFC_layer_name);
				std::vector<float> LPFCLayerValues;

				if (vector_of_layers.size() > 0)
				{
					auto layer = vector_of_layers[0];
					all_points.GetLayerVals(LPFCLayerValues, *layer);

					// create data layers
					auto RegionsLayer = cloud.CreateLayer(L"Regions", 0.0);
					std::vector<float> RegionsIDs;

					auto SegmentsLayer = cloud.CreateLayer(L"Segments", 0.0);
					std::vector<float> SegmentsIDs;

					auto senttostackLayer = cloud.CreateLayer(L"Senttostack", 0.0);
					std::vector<float> val_senttostack(all_points.size(), 0.0);

					std::vector<Math::Vector3D> regions_normals;
					std::vector<float> regions_sizes;

					// Pick random points from the cloud as initial seeds
					Data::Clouds::RandomSearchKernel RandomSearch(number_of_points, false);
					Data::Clouds::PointsRange random_points_range;
					cloud.GetAccess().FindPoints(RandomSearch, random_points_range);

					std::vector<Data::Clouds::Point3D> random_points;
					random_points_range.GetXYZ(random_points);

					std::vector<Data::Clouds::Point3D> rp_normals;
					random_points_range.GetNormals(rp_normals);

					std::stack<Math::Point3D> seed_points;

					float new_region_id = 1.0;

					// Creating regions
					for (auto i = 0; i < random_points_range.size(); i++)
					{
						std::vector<float> regions_ids;
						random_points_range.GetLayerVals(regions_ids, *RegionsLayer);

						std::vector<float> rp_LPFC_values;
						random_points_range.GetLayerVals(rp_LPFC_values, *layer);

						if (regions_ids[i] > 0 || rp_LPFC_values[i] > min_LPFC) // if point already in a region or its LPFC parameter value is too large -> continue
							continue;

						auto initial_seed = random_points[i].cast<Real>();
						seed_points.push(initial_seed);

						regions_ids[i] = new_region_id; // set region id of the initial seed point here
						random_points_range.SetLayerVals(regions_ids, *RegionsLayer);

						Math::Vector3D region_normal_vect = rp_normals[i].cast<Real>(); // set seed point's normal as region's normal

						std::vector<Math::Point3D> cluster; // store points added to the region
						std::vector<Math::Vector3D> region_points_normals; // store the normal vectors of the points added to the region

						float saved_cluster_size = 0;
						bool update_region_normal = true;

						// Region growing
						while (!seed_points.empty())
						{
							Math::Point3D seed = seed_points.top();
							seed_points.pop();

							// create local neighbourhood around seed point using KNN
							Data::Clouds::KNNSearchKernel search(seed, KNN_size);
							Data::Clouds::PointsRange neighbourhood;
							cloud.GetAccess().FindPoints(search, neighbourhood);

							std::vector<Data::Clouds::Point3D> neighbourhood_points;
							neighbourhood.GetXYZ(neighbourhood_points);

							std::vector<Data::Clouds::Point3D> neighbourhood_normals;
							neighbourhood.GetNormals(neighbourhood_normals);

							std::vector<float> nb_region_ids;
							std::vector<float> nb_senttostack;

							neighbourhood.GetLayerVals(nb_region_ids, *RegionsLayer);
							neighbourhood.GetLayerVals(nb_senttostack, *senttostackLayer);

							for (auto ii = 0; ii < neighbourhood.size(); ii++)
							{
								if (nb_region_ids[ii] > 0)
									continue;

								Math::Point3D nb_point = neighbourhood_points[ii].cast<Real>();
								Math::Vector3D nb_point_normal_vect = neighbourhood_normals[ii].cast<Real>();

								Real angle = Math::CalcAngleBetweenTwoVectors(region_normal_vect, nb_point_normal_vect) * 180 / PI; // angle in degrees

								if (angle <= threshold_angle)
								{
									// add point to the region
									nb_region_ids[ii] = new_region_id; // set point's region id

									cluster.push_back(nb_point); // add point to cluster
									region_points_normals.push_back(nb_point_normal_vect);

									if (nb_senttostack[ii] == 0.0 && nb_point != initial_seed)
									{
										seed_points.push(nb_point);
										nb_senttostack[ii] = 1;
									}

									if (!update_region_normal)
										continue;

									float growth_factor = 100 * (cluster.size() - saved_cluster_size) / saved_cluster_size; // current region growth factor

									if ((cluster.size() >= 3) && (growth_factor >= threshold_growth_factor)) // at least 3 points required for plane fitting
									{
										Math::Vector3D new_region_normal;

										// calculate new region normal by fitting a plane to the points in the current region
										auto plane = Math::CalcBestPlane3D(cluster.begin(), cluster.end());
										Math::Vector3D surface_normal = plane.normal();

										bool same_direction = is_same_direction(surface_normal, region_normal_vect); // check if surface normal and region normal are in the same direction
										if (same_direction)
											new_region_normal = surface_normal;
										else
											new_region_normal = (-1) * surface_normal; // flip vector

										bool new_region_normal_is_ok = true;
										for (auto region_point_normal : region_points_normals) // check if new region normal is ok - every normal vector in the region must be within the threshold angle
										{
											Real angle = Math::CalcAngleBetweenTwoVectors(region_point_normal, new_region_normal) * 180 / PI; // angle in degrees
											if (angle > threshold_angle)
											{
												new_region_normal_is_ok = false;
												break;
											}
										}

										if (new_region_normal_is_ok)
										{
											region_normal_vect = new_region_normal;
											saved_cluster_size = cluster.size();
										}
										else
											update_region_normal = false;
									}
								}
							}

							neighbourhood.SetLayerVals(nb_region_ids, *RegionsLayer);
							neighbourhood.SetLayerVals(nb_senttostack, *senttostackLayer);
						}

						float region_size = cluster.size();

						//OGX_LINE.Format(Info, L"Region size: %f", region_size);

						regions_sizes.push_back(region_size);
						regions_normals.push_back(region_normal_vect);

						new_region_id++;
					}

					all_points.GetLayerVals(RegionsIDs, *RegionsLayer);
					all_points.GetLayerVals(SegmentsIDs, *SegmentsLayer);

					// Region grouping - segmentation
					std::vector<std::pair<int, int>> regions; // first - region id, second - region size

					for (auto i = 0; i < regions_sizes.size(); i++)
						regions.push_back(std::make_pair(i + 1, regions_sizes[i]));

					// sorting regions by size
					sort(regions.begin(), regions.end(), sortbysecdesc); // regions sorted by their size

					// segmentation
					std::vector<Math::Vector3D> segments_normals;

					for (auto region : regions)
					{
						int current_region_id = region.first; // i + 1
						int current_region_size = region.second;

						//OGX_LINE.Format(Info, L"segment_id: %f", current_region_id);

						Math::Vector3D current_region_normal = regions_normals[current_region_id - 1];

						// define threshold angle by size
						Real min_angle = 1;
						if (current_region_size <= SR_size_threshold)
							min_angle = SR_threshold_angle;

						if (min_angle <= 0)
							min_angle = 1;

						// search for the best segment
						bool segment_found = false;
						int segment_id = 0.0;

						for (auto ii = 0; ii < segments_normals.size(); ii++)
						{
							Math::Vector3D segment_normal = segments_normals[ii];

							Real angle = Math::CalcAngleBetweenTwoVectors(current_region_normal, segment_normal) * 180 / PI; // angle in degrees

							if (angle < min_angle)
							{
								min_angle = angle;
								segment_id = ii + 1;
								segment_found = true;
							}
						}

						if (segment_found)
						{
							// add region to segment
							for (auto ii = 0; ii < all_points.size(); ii++)
								if (RegionsIDs[ii] == current_region_id)
									SegmentsIDs[ii] = segment_id;
						}
						else
						{
							// create new segment
							segment_id = segments_normals.size() + 1;

							for (auto ii = 0; ii < all_points.size(); ii++)
								if (RegionsIDs[ii] == current_region_id)
									SegmentsIDs[ii] = segment_id;

							segments_normals.push_back(current_region_normal);
						}
					}

					all_points.SetLayerVals(SegmentsIDs, *SegmentsLayer);

					// Saving viewing vectors to data layers
					auto ViewingVectorsXLayer = cloud.CreateLayer(L"View_X", 0.0); // stores the x component of the viewing vector
					std::vector<float> val_X(all_points.size(), 0.0);

					auto ViewingVectorsYLayer = cloud.CreateLayer(L"View_Y", 0.0); // stores the y component of the viewing vector
					std::vector<float> val_Y(all_points.size(), 0.0);

					auto ViewingVectorsZLayer = cloud.CreateLayer(L"View_Z", 0.0); // stores the z component of the viewing vector
					std::vector<float> val_Z(all_points.size(), 0.0);

					auto ClassifiedLayer = cloud.CreateLayer(L"Classified", 0.0); // stores info whether point has been classified to a segment
					std::vector<float> ClassifiedValues(all_points.size(), 0.0); // 1 - classified, 0 - not classified

					float count = 0.0;

					for (auto i = 0; i < all_points.size(); i++)
					{
						if (SegmentsIDs[i] == 0.0)
						{
							count++;
							continue;
						}

						ClassifiedValues[i] = 1; // point classified as segmented

						int vector_id = SegmentsIDs[i] - 1;

						val_X[i] = (-1) * segments_normals[vector_id].x();
						val_Y[i] = (-1) * segments_normals[vector_id].y();
						val_Z[i] = (-1) * segments_normals[vector_id].z();
					}

					all_points.SetLayerVals(val_X, *ViewingVectorsXLayer);
					all_points.SetLayerVals(val_Y, *ViewingVectorsYLayer);
					all_points.SetLayerVals(val_Z, *ViewingVectorsZLayer);

					all_points.SetLayerVals(ClassifiedValues, *ClassifiedLayer);

					float number_of_segments = segments_normals.size();
					OGX_LINE.Format(Info, L"Number of segments: %f", number_of_segments);

					float segmentation_quality = 100 * (all_points.size() - count) / all_points.size();
					OGX_LINE.Format(Info, L"Number of unclassified points: %f, Segmentation quality: %f", count, segmentation_quality);

					end = std::chrono::steady_clock::now();
					float elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
					OGX_LINE.Format(Info, L"Total elapsed time: %f ms", elapsed_time);
				}
				else
					ReportError(L"There are no layers");

			}, thread_count); // run with given number of threads, optional parameter, if not given will run in current threadlefrw0o
	}
};

/** \struct
\brief Normal vector histogram segmentation method

\Divides point cloud into segments of similarly oriented surfaces using a normal vector histogram

*/

struct SegmentationByNormalVectorHistogram : public ogx::Plugin::EasyMethod
{
	// parameters
	Data::ResourceID m_node_id;
	Real threshold_angle;
	Real PPU;

	// constructor
	SegmentationByNormalVectorHistogram() : EasyMethod(L"Micha³ Kossakowski", L"Segmentation using a normal vector histogram.")
	{
	}

	// add input/output parameters
	virtual void DefineParameters(ParameterBank& bank)
	{
		bank.Add(L"node_id", m_node_id, L"ID of node containing point cloud(s)").AsNode();
		bank.Add(L"threshold_angle", threshold_angle = 15, L"Acceptable deviation of surface orientation [degrees]"); // angle in degrees
		bank.Add(L"PPU", PPU = 100, L"Percentage of Points Used - percentage of total number of points contained in the most popular histogram peaks [%]");

	}

	virtual void Run(Context& context)
	{
		auto subtree = context.Project().TransTreeFindNode(m_node_id);
		// report error if given node was not found, this will stop execution of algorithm
		if (!subtree) ReportError(L"Node not found");

		// run with number of threads available on current machine, optional
		auto const thread_count = std::thread::hardware_concurrency();

		// perform calculations for each cloud in given subtree
		Clouds::ForEachCloud(*subtree, [&](Data::Clouds::ICloud& cloud, Data::Nodes::ITransTreeNode& node)
			{
				// start timer
				auto start = std::chrono::steady_clock::now();

				//access points in the cloud
				Data::Clouds::PointsRange all_points;
				cloud.GetAccess().GetAllPoints(all_points);

				std::vector<Data::Clouds::Point3D> points_xyz;
				all_points.GetXYZ(points_xyz);

				std::vector<Data::Clouds::Point3D> normals;
				all_points.GetNormals(normals);

				// create data layers
				auto BucketsLayer = cloud.CreateLayer(L"Buckets", 0.0);
				std::vector<float> BucketsIDs;

				auto FrequencyLayer = cloud.CreateLayer(L"Frequency", 0.0);
				std::vector<float> FrequencyValues;

				auto SegmentsLayer = cloud.CreateLayer(L"Segments", 0.0);
				std::vector<float> SegmentsIDs(all_points.size(), 0.0);

				auto PhiLayer = cloud.CreateLayer(L"Phi", 0.0);
				std::vector<float> PhiValues(all_points.size(), 0.0);

				auto ThetaLayer = cloud.CreateLayer(L"Theta", 0.0);
				std::vector<float> ThetaValues(all_points.size(), 0.0);

				float histogram_bucket_size = 2 * threshold_angle;

				// calculate number of buckets in the phi and theta axes: theta - n, phi - m
				const int n = std::ceil((1.0) * 360 / histogram_bucket_size);
				const int m = std::ceil((1.0) * 180 / histogram_bucket_size);

				histogram_bucket_size = (1.0) * 180 / m;

				// create histogram
				std::map<int, int> histogram; // first - bucket_id, second - frequency

				for (auto i = 0; i < all_points.size(); i++)
				{
					auto normal = normals[i];

					float bucket_id = 0.0;

					float r = 0.0;
					float phi = 0.0;
					float theta = 0.0;

					// normal vector components
					Real x = normal.x();
					Real y = normal.y();
					Real z = normal.z();

					// calculating spherical coordinates
					r = sqrt(x * x + y * y + z * z);
					phi = acos(z / r) * 180 / PI; // (0; 180)
					theta = atan2(y, x) * 180 / PI; // (-90; 90)

					if (theta < 0)
						theta = 360 + theta;

					// assigning normal vector to the appropriate bucket
					int phi_id = phi / histogram_bucket_size + 1;
					int theta_id = theta / histogram_bucket_size + 1;

					if (phi == 180)
						phi_id = m;
					if (theta == 360)
						theta_id = n;

					bucket_id = n * (phi_id - 1) + theta_id;

					// update frequency value in histogram
					histogram[bucket_id]++;

					BucketsIDs.push_back(bucket_id);

					PhiValues[i] = phi;
					ThetaValues[i] = theta;
				}

				all_points.SetLayerVals(BucketsIDs, *BucketsLayer);

				all_points.SetLayerVals(PhiValues, *PhiLayer);
				all_points.SetLayerVals(ThetaValues, *ThetaLayer);

				// set frequency layer values
				for (auto i = 0; i < all_points.size(); i++)
				{
					float frequency = histogram[BucketsIDs[i]];
					FrequencyValues.push_back(frequency);
				}

				all_points.SetLayerVals(FrequencyValues, *FrequencyLayer);

				// sort the histogram by the frequency value in descending order
				std::vector<std::pair<int, int>> sorted_histogram; // first - bucket_id, second - frequency

				for (auto& it : histogram)
					sorted_histogram.push_back(it);

				sort(sorted_histogram.begin(), sorted_histogram.end(), sortbysecdesc);

				// Viewing vectors calculation and segmentation
				int number_of_points_covered = 0;
				int number_of_points_limit = (PPU / 100) * all_points.size();

				std::vector<Math::Vector3D> measuring_directions;

				int segment_id = 1;

				for (auto& it : sorted_histogram)
				{
					if (number_of_points_covered >= number_of_points_limit)
						continue;

					int bucket_id = it.first;
					int frequency = it.second;

					number_of_points_covered += frequency;

					std::vector<Math::Point3D> cluster;

					Math::Vector3D bucket_orientation_vector; // defines the direction of the normals in current bucket
					Math::Vector3D measuring_direction; // defines the direction of measurement for points in curtrent bucket

					// Assigning points to segment
					for (auto i = 0; i < all_points.size(); i++)
						if (BucketsIDs[i] == bucket_id)
						{
							SegmentsIDs[i] = segment_id;

							if (cluster.size() == 0) // for the first point save its normal vector
								bucket_orientation_vector = normals[i].cast<Real>();

							Math::Point3D point = points_xyz[i].cast<Real>();
							cluster.push_back(point);
						}

					// Calculate measuring direction by plane fitting
					auto plane = Math::CalcBestPlane3D(cluster.begin(), cluster.end());
					Math::Vector3D surface_normal = plane.normal();

					Real angle = Math::CalcAngleBetweenTwoVectors(surface_normal, bucket_orientation_vector);
					if (angle > 90)
						measuring_direction = surface_normal; // flip the calculated vector if necessary
					else
						measuring_direction = (-1) * surface_normal;

					measuring_directions.push_back(measuring_direction);

					segment_id++;
				}

				all_points.SetLayerVals(SegmentsIDs, *SegmentsLayer);

				float count = 0.0;

				// Saving normalized viewing vectors to data layers
				auto ViewingVectorsXLayer = cloud.CreateLayer(L"View_X", 0.0); // stores the x component of the viewing vector
				std::vector<float> val_X(all_points.size(), 0.0);

				auto ViewingVectorsYLayer = cloud.CreateLayer(L"View_Y", 0.0); // stores the y component of the viewing vector
				std::vector<float> val_Y(all_points.size(), 0.0);

				auto ViewingVectorsZLayer = cloud.CreateLayer(L"View_Z", 0.0); // stores the z component of the viewing vector
				std::vector<float> val_Z(all_points.size(), 0.0);

				auto ClassifiedLayer = cloud.CreateLayer(L"Classified", 0.0); // stores info whether point has been classified to a segment
				std::vector<float> ClassifiedValues(all_points.size(), 0.0); // 1 - classified, 0 - not classified

				for (auto i = 0; i < all_points.size(); i++)
				{
					if (SegmentsIDs[i] == 0.0)
					{
						count++;
						continue;
					}

					ClassifiedValues[i] = 1;

					// saving viewing vector
					int vector_id = SegmentsIDs[i] - 1;

					// saving bucket viewing vector components
					val_X[i] = measuring_directions[vector_id].x(); // X component
					val_Y[i] = measuring_directions[vector_id].y(); // Y component
					val_Z[i] = measuring_directions[vector_id].z(); // Z component
				}

				all_points.SetLayerVals(val_X, *ViewingVectorsXLayer);
				all_points.SetLayerVals(val_Y, *ViewingVectorsYLayer);
				all_points.SetLayerVals(val_Z, *ViewingVectorsZLayer);

				all_points.SetLayerVals(ClassifiedValues, *ClassifiedLayer);

				float segmentation_quality = 100 * (all_points.size() - count) / all_points.size();

				OGX_LINE.Format(Info, L"Number of unclassified points: %f, Segmentation quality: %f", count, segmentation_quality);

				float number_of_segments = 1.0 * segment_id - 1;

				OGX_LINE.Format(Info, L"Number of segments: %f", number_of_segments);

				// stop timer
				auto end = std::chrono::steady_clock::now();
				float elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

				OGX_LINE.Format(Info, L"Elapsed time: %f ms", elapsed_time);

			}, thread_count); // run with given number of threads, optional parameter, if not given will run in current threadlefrw0o
	}
};

// Export methods

OGX_EXPORT_METHOD(Local_plane_fitting_coefficient)
OGX_EXPORT_METHOD(SegmentationByRegionGrowingWithLPFC)
OGX_EXPORT_METHOD(SegmentationByNormalVectorHistogram)