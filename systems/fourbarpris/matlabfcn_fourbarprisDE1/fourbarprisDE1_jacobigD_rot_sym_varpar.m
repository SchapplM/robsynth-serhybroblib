% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% fourbarprisDE1
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% qJD [1x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[GK,GP,HP]';
% 
% Output:
% JgD_rot [3x1]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 09:10
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = fourbarprisDE1_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),uint8(0),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbarprisDE1_jacobigD_rot_sym_varpar: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbarprisDE1_jacobigD_rot_sym_varpar: qJD has to be [1x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'fourbarprisDE1_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisDE1_jacobigD_rot_sym_varpar: pkin has to be [3x1] (double)');
JgD_rot=NaN(3,1);
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 09:10:20
	% EndTime: 2020-05-07 09:10:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0; 0; 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 09:10:20
	% EndTime: 2020-05-07 09:10:20
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (92->17), mult. (47->16), div. (3->3), fcn. (6->2), ass. (0->11)
	t64 = qJ(1) + pkin(3);
	t69 = -pkin(2) + t64;
	t61 = pkin(1) + t69;
	t70 = -pkin(2) - t64;
	t62 = pkin(1) + t70;
	t71 = t61 * t62;
	t60 = pkin(1) - t69;
	t68 = t60 * t71;
	t59 = pkin(1) - t70;
	t67 = t59 * t68;
	t1 = [0; 0; (-0.2e1 + (-0.1e1 / t64 ^ 2 + (-t68 + (t71 + (t61 - t62) * t60) * t59) / t67 / t64 / 0.2e1) * (pkin(1) ^ 2 - pkin(2) ^ 2 - qJ(1) ^ 2 + (-2 * qJ(1) - pkin(3)) * pkin(3))) * (-t67) ^ (-0.1e1 / 0.2e1) * qJD(1);];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 09:10:20
	% EndTime: 2020-05-07 09:10:20
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (92->17), mult. (47->16), div. (3->3), fcn. (6->2), ass. (0->11)
	t46 = qJ(1) + pkin(3);
	t51 = -pkin(2) + t46;
	t43 = pkin(1) + t51;
	t52 = -pkin(2) - t46;
	t44 = pkin(1) + t52;
	t53 = t43 * t44;
	t42 = pkin(1) - t51;
	t50 = t42 * t53;
	t41 = pkin(1) - t52;
	t49 = t41 * t50;
	t1 = [0; 0; (-0.2e1 + (-0.1e1 / t46 ^ 2 + (-t50 + (t53 + (t43 - t44) * t42) * t41) / t49 / t46 / 0.2e1) * (pkin(1) ^ 2 - pkin(2) ^ 2 - qJ(1) ^ 2 + (-2 * qJ(1) - pkin(3)) * pkin(3))) * (-t49) ^ (-0.1e1 / 0.2e1) * qJD(1);];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 09:10:20
	% EndTime: 2020-05-07 09:10:20
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (67->12), mult. (24->9), div. (0->1), fcn. (4->2), ass. (0->11)
	t56 = (pkin(3) + qJ(1));
	t54 = -pkin(2) + t56;
	t48 = (pkin(1) + t54);
	t55 = -pkin(2) - t56;
	t49 = (pkin(1) + t55);
	t57 = (t48 * t49);
	t47 = (pkin(1) - t54);
	t53 = (t47 * t57);
	t46 = (pkin(1) - t55);
	t52 = (t46 * t53);
	t1 = [0; 0; (-2 - t56 * (-t53 + (t57 + (t48 - t49) * t47) * t46) / t52) * (-t52) ^ (-0.1e1 / 0.2e1) * qJD(1);];
	JgD_rot = t1;
end