% Jacobian time derivative of explicit kinematic constraints of
% fourbarprisIC
% with respect to passive joint coordinates
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[GK,GP,HP]';
% 
% Output:
% PhiD_p [(no of constraints)x(no. of passive joints)] 

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 09:59
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function PhiD_p = fourbarprisIC_kinconstr_impl_pas_jacobianD_mdh_sym_varpar(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbarprisIC_kinconstr_impl_pas_jacobianD_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'fourbarprisIC_kinconstr_impl_pas_jacobianD_mdh_sym_varpar: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisIC_kinconstr_impl_pas_jacobianD_mdh_sym_varpar: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From kinconstr_impl_passive_jacobianD_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:59:30
% EndTime: 2020-05-07 09:59:30
% DurationCPUTime: 0.01s
% Computational Cost: add. (6->5), mult. (10->8), div. (0->0), fcn. (6->4), ass. (0->5)
t9 = pkin(2) * qJD(3);
t8 = qJD(1) * (pkin(3) + qJ(2));
t7 = cos(qJ(1));
t6 = sin(qJ(1));
t1 = [t6 * qJD(2) + t7 * t8, -cos(qJ(3)) * t9, 0; -qJD(2) * t7 + t6 * t8, -sin(qJ(3)) * t9, 0; 0, 0, 0;];
PhiD_p = t1;
